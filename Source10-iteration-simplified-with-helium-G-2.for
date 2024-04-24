      SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd,
     1 rpl, ddsddt, drplde, drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c  WRITE (6,*) '
c  NOTE:  MODIFICATIONS TO *UMAT FOR ABAQUS VERSION 5.3 (14 APR '94)
c 
c  (1)  The list of variables above defining the *UMAT subroutine, 
c  and the first (standard) block of variables dimensioned below, 
c  have variable names added compared to earlier ABAQUS versions. 
c
c  (2)  The statement: include 'aba_param.inc' must be added as below.
c
c  (3)  As of version 5.3, ABAQUS files use double precision only.
c  The file aba_param.inc has a line "implicit real*8" and, since 
c  it is included in the main subroutine, it will define the variables
c  there as double precision.  But other subroutines still need the
c  definition "implicit real*8" since there may be variables that are
c  not passed to them through the list or common block.
c
c  (4)  This is current as of version 5.6 of ABAQUS.
c
c  (5)  Note added by J. W. Kysar (4 November 1997).  This UMAT has been
c   modified to keep track of the cumulative shear strain in each 
c   individual slip system.  This information is needed to correct an
c   error in the implementation of the Bassani and Wu hardening law.
c   Any line of code which has been added or modified is preceded
c   immediately by a line beginning CFIXA and succeeded by a line 
c   beginning CFIXB.  Any comment line added or modified will begin 
c   with CFIX.
c
c   The hardening law by Bassani and Wu was implemented incorrectly.
c   This law is a function of both hyperbolic secant squared and hyperbolic
c   tangent.  However, the arguments of sech and tanh are related to the *total*
c   slip on individual slip systems.  Formerly, the UMAT implemented this
c   hardening law by using the *current* slip on each slip system.  Therein 
c   lay the problem. The UMAT did not restrict the current slip to be a 
c   positive value.  So when a slip with a negative sign was encountered, the 
c   term containing tanh led to a negative hardening rate (since tanh is an 
c   odd function).

c   The UMAT has been fixed by adding state variables to keep track of the 
c   *total* slip on each slip system by integrating up the absolute value 
c   of slip rates for each individual slip system.  These "solution dependent 
c   variables" are available for postprocessing.  The only required change 
c   in the input file is that the DEPVAR command must be changed.
c 
C-----  Use single precision on Cray by
C     (1) deleting the statement "IMPLICIT*8 (A-H,O-Z)";
C     (2) changing "REAL*8 FUNCTION" to "FUNCTION";
C     (3) changing double precision functions DSIGN to SIGN.
C
C-----  Subroutines:
C
C       ROTATION     -- forming rotation matrix, i.e. the direction 
C                       cosines of cubic crystal [100], [010] and [001]
C                       directions in global system at the initial 
C                       state
C
C       SLIPSYS      -- calculating number of slip systems, unit 
C                       vectors in slip directions and unit normals to 
C                       slip planes in a cubic crystal at the initial 
C                       state
C
C       GSLPINIT     -- calculating initial value of current strengths 
C                       at initial state
C
C       STRAINRATE   -- based on current values of resolved shear 
C                       stresses and current strength, calculating 
C                       shear strain-rates in slip systems
C
C       LUDCMP       -- LU decomposition
C
C       LUBKSB       -- linear equation solver based on LU 
C                       decomposition method (must call LUDCMP first)


C-----  Function subprogram:

C       F -- shear strain-rates in slip systems


C-----  Variables:
C
C       STRESS -- stresses (INPUT & OUTPUT)
C                 Cauchy stresses for finite deformation
C       STATEV -- solution dependent state variables (INPUT & OUTPUT)
C       DDSDDE -- Jacobian matrix (OUTPUT)

C-----  Variables passed in for information:
C
C       STRAN  -- strains
C                 logarithmic strain for finite deformation 
C                 (actually, integral of the symmetric part of velocity
C                  gradient with respect to time)
C       DSTRAN -- increments of strains
C       CMNAME -- name given in the *MATERIAL option
C       NDI    -- number of direct stress components
C       NSHR   -- number of engineering shear stress components
C       NTENS  -- NDI+NSHR
C       NSTATV -- number of solution dependent state variables (as 
C                 defined in the *DEPVAR option)
C       PROPS  -- material constants entered in the *USER MATERIAL 
C                 option
C       NPROPS -- number of material constants
C

C-----  This subroutine provides the plastic constitutive relation of 
C     single crystals for finite element code ABAQUS.  The plastic slip
C     of single crystal obeys the Schmid law.  The program gives the 
C     choice of small deformation theory and theory of finite rotation 
C     and finite strain.
C       The strain increment is composed of elastic part and plastic 
C     part.  The elastic strain increment corresponds to lattice 
C     stretching, the plastic part is the sum over all slip systems of 
C     plastic slip.  The shear strain increment for each slip system is
C     assumed a function of the ratio of corresponding resolved shear 
C     stress over current strength, and of the time step.  The resolved
C     shear stress is the double product of stress tensor with the slip
C     deformation tensor (Schmid factor), and the increment of current 
C     strength is related to shear strain increments over all slip 
C     systems through self- and latent-hardening functions.

C-----  The implicit integration method proposed by Peirce, Shih and 
C     Needleman (1984) is used here.  The subroutine provides an option
C     of iteration to solve stresses and solution dependent state 
C     variables within each increment.

C-----  The present program is for a single CUBIC crystal.  However, 
C     this code can be generalized for other crystals (e.g. HCP, 
C     Tetragonal, Orthotropic, etc.).  Only subroutines ROTATION and 
C     SLIPSYS need to be modified to include the effect of crystal 
C     aspect ratio.
C

C-----  Important notice:
C
C     (1) The number of state variables NSTATV must be larger than (or 
CFIX      equal to) TEN (10) times the total number of slip systems in
C         all sets, NSLPTL, plus FIVE (5)
CFIX           NSTATV >= 10 * NSLPTL + 5
C         Denote s as a slip direction and m as normal to a slip plane.
C         Here (s,-m), (-s,m) and (-s,-m) are NOT considered 
C         independent of (s,m).  The number of slip systems in each set
C         could be either 6, 12, 24 or 48 for a cubic crystal, e.g. 12 
C         for {110}<111>.
C
C         Users who need more parameters to characterize the 
C         constitutive law of single crystal, e.g. the framework 
C         proposed by Zarka, should make NSTATV larger than (or equal 
C         to) the number of those parameters NPARMT plus nine times 
C         the total number of slip systems, NSLPTL, plus five
CFIX           NSTATV >= NPARMT + 10 * NSLPTL + 5
C
C     (2) The tangent stiffness matrix in general is not symmetric if 
C         latent hardening is considered.  Users must declare "UNSYMM" 
C         in the input file, at the *USER MATERIAL card.
C

      PARAMETER (ND=15)
C-----  The parameter ND determines the dimensions of the arrays in 
C     this subroutine.  The current choice 150 is a upper bound for a 
C     cubic crystal with up to three sets of slip systems activated.  
C     Users may reduce the parameter ND to any number as long as larger
C     than or equal to the total number of slip systems in all sets.  
C     For example, if {110}<111> is the only set of slip system 
C     potentially activated, ND could be taken as twelve (12).  
c
      include 'aba_param.inc'

c
      CHARACTER*8 CMNAME
      CHARACTER*80 PRTSTR

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), SLPDIR(3,ND),
     2          SLPNOR(3,ND), SLPDEF(6,ND), SLPDEFNOR(6,ND),
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND), 
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3), 
     5          FSLIP(ND), DFDXTAU(ND), DFDXRHO(ND,ND), DFDXSIG(ND),
     6          DDEMSD(6,ND), DDEMSDNOR(6,ND), DDGDDE(ND,6), 
     7          DSTRES(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DRHOSP(ND), DKSP(ND),
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3),
     1          ELASTRAN(6), PLASTRAN(6), DELATS(6), DPLATS(6),
     2          ROTSTIFF(16), DPLATSM(6)


      DIMENSION DDRDDGA(ND,ND), DDRDDGC(ND,ND), DRHORK45(ND,4),
     2          SLIPDM(ND,ND), SLIPGM(ND,ND), SLIPHM(ND,ND)


      DIMENSION FSLIP1(ND), STRES1(6),  
     2          SPNOR1(3,ND), SPDIR1(3,ND), DDSDE1(6,6),
     3          DSOLD(6), DSPNRO(3,ND), DSPDRO(3,ND), GAMMA1(ND),
     4          TAUSP1(ND), RHOSP1(ND), GSLP1(ND), SIGSP1(ND),
     5          TAUSPEFF1(ND), SIGSPEFF1(ND),
     6          DGAMOD(ND), DTAUOD(ND), DRHOOD(ND), SGAMMA1(ND),
     7          ELASTRAN1(6),PLASTRAN1(6), DELATSOD(6), DPLATSOD(6),
     8          TOTALSTRAN1(6)

      DIMENSION TAUSPEFF(ND), DEVSIG(3,3), SLPDEFEFF33(3,3,ND), 
     2          SLPDEFEFF(6,ND), SLPSPNEFF(3,ND), DFDXTAUEFF(ND), 
     3          DTAUEFFDTAU(ND), DTAUEFFDFVOID(ND), DFDXSIGEFF(ND),
     4          DPHIDTAU(ND), DPHIDTAUEFF(ND), DPHIDFVOIDEFF(ND), 
     5          DDFVOIDDDG(ND), DDEMSDEFF(6,ND), SIGSPEFF(ND),
     6          IDTAUEFF(ND), DDFVOIDDDE(6), SLPDEFNOREFF33(3,3,ND),
     7          SLPDEFNOREFF(6,ND), DDEMSDNOREFF(6,ND), 
     8          DDEMSDEFF0(6,ND)
     
      REAL*8 Weibull(50000) 
     
C-----  NSLIP  -- number of slip systems in each set
C-----  SLPDIR -- slip directions (unit vectors in the initial state)
C-----  SLPNOR -- normals to slip planes (unit normals in the initial 
C                 state)
C-----  SLPDEF -- slip deformation tensors (Schmid factors)
C                 SLPDEF(1,i) -- SLPDIR(1,i)*SLPNOR(1,i)
C                 SLPDEF(2,i) -- SLPDIR(2,i)*SLPNOR(2,i)
C                 SLPDEF(3,i) -- SLPDIR(3,i)*SLPNOR(3,i)
C                 SLPDEF(4,i) -- SLPDIR(1,i)*SLPNOR(2,i)+
C                                SLPDIR(2,i)*SLPNOR(1,i)
C                 SLPDEF(5,i) -- SLPDIR(1,i)*SLPNOR(3,i)+
C                                SLPDIR(3,i)*SLPNOR(1,i)
C                 SLPDEF(6,i) -- SLPDIR(2,i)*SLPNOR(3,i)+
C                                SLPDIR(3,i)*SLPNOR(2,i)
C                 where index i corresponds to the ith slip system
C-----  SLPSPN -- slip spin tensors (only needed for finite rotation)
C                 SLPSPN(1,i) -- [SLPDIR(1,i)*SLPNOR(2,i)-
C                                 SLPDIR(2,i)*SLPNOR(1,i)]/2
C                 SLPSPN(2,i) -- [SLPDIR(3,i)*SLPNOR(1,i)-
C                                 SLPDIR(1,i)*SLPNOR(3,i)]/2
C                 SLPSPN(3,i) -- [SLPDIR(2,i)*SLPNOR(3,i)-
C                                 SLPDIR(3,i)*SLPNOR(2,i)]/2
C                 where index i corresponds to the ith slip system
C-----  DSPDIR -- increments of slip directions
C-----  DSPNOR -- increments of normals to slip planes
C
C-----  DLOCAL -- elastic matrix in local cubic crystal system
C-----  D      -- elastic matrix in global system
C-----  ROTD   -- rotation matrix transforming DLOCAL to D
C
C-----  ROTATE -- rotation matrix, direction cosines of [100], [010] 
C                 and [001] of cubic crystal in global system
C
C-----  FSLIP  -- shear strain-rates in slip systems
C----------------------------------------------------------------------
C-----  DFDXTAU -- derivatives of FSLIP w.r.t x=TAUSLP, where 
C                 TAUSLP is the resolved shear stress
C-----  DFDXRHO -- derivatives of FSLIP w.r.t x=RHOSLP, where 
C                 RHOSLP is the dislocation density
C
C-----  DDEMSD -- double dot product of the elastic moduli tensor with 
C                 the slip deformation tensor plus, only for finite 
C                 rotation, the dot product of slip spin tensor with 
C                 the stress
C-----  DDGDDE -- derivatice of the shear strain increments in slip 
C                 systems w.r.t. the increment of strains
C-----  DDRDDG -- derivatice of the dislocation increments in slip 
C                 systems w.r.t. the increment of slip value
C
C-----  DSTRES -- Jaumann increments of stresses, i.e. corotational 
C                 stress-increments formed on axes spinning with the 
C                 material
C-----  DELATS -- strain-increments associated with lattice stretching
C                 DELATS(1) - DELATS(3) -- normal strain increments
C                 DELATS(4) - DELATS(6) -- engineering shear strain 
C                                          increments
C-----  DSPIN  -- spin-increments associated with the material element
C                 DSPIN(1) -- component 12 of the spin tensor
C                 DSPIN(2) -- component 31 of the spin tensor
C                 DSPIN(3) -- component 23 of the spin tensor
C
C-----  DVGRAD -- increments of deformation gradient in the current 
C                 state, i.e. velocity gradient times the increment of 
C                 time
C
C-----  DGAMMA -- increment of shear strains in slip systems
C-----  DTAUSP -- increment of resolved shear stresses in slip systems 
C-----  DRHOSP -- increment of dislocation density in slip systems
C-----  DKSP   -- increment of K in slip systems
C
C-----  Arrays for iteration:
C
C            FSLIP1, STRES1, DDSDE1, DSOLD ,
C-----       DGAMOD, DTAUOD, DRHOOD, DKOD, DSPNRO, DSPDRO 
C            GAMMA1, TAUSP1, RHOSP1, KSP1, SPNOR1, SPDIR1, GSLP1
C
C
C-----  Solution dependent state variable STATEV:
C            Denote the number of total slip systems by NSLPTL, which 
C            will be calculated in this code.
C
C       Array STATEV:
C       1          - NSLPTL    :  current strength in slip systems
C       NSLPTL+1   - 2*NSLPTL  :  shear strain in slip systems
C       2*NSLPTL+1 - 3*NSLPTL  :  resolved shear stress in slip systems
C
C       3*NSLPTL+1 - 6*NSLPTL  :  current components of normals to slip
C                                 slip planes
C       6*NSLPTL+1 - 9*NSLPTL  :  current components of slip directions
C
CFIX    9*NSLPTL+1 - 10*NSLPTL :  total cumulative shear strain on each 
CFIX                              slip system (sum of the absolute 
CFIX                              values of shear strains in each slip 
CFIX                              system individually)
CFIX
CFIX    10*NSLPTL+1 - 11*NSLPTL : resolved normal stress in slip sysytems
CFIX
C
CFIX    11*NSLPTL+1-12*NSLPTL+6 : elastic strain components

CFIX    11*NSLPTL+7-12*NSLPTL+12 : plastic strain components

CFIX    11*NSLPTL+13-12*NSLPTL+18 : total strain components
C
C
C
CFIX    11*NSLPTL+19             :  equivalent plastic strain
C


CFIX    11*NSLPTL+20             : total dislocation density on all slip
CFIX                               systems
C
CFIX    11*NSLPTL+21-11*NSLPTL+26: crystal orientation after rigid body rotation



C       11*NSLPTL+27             : equivalent strain


CFIX    11*NSLPTL+28 - 12*NSLPTL+27 : dislocation density in each slip system
CFIX                                  individually
CFIX    12*NSLPTL+28 - 13*NSLPTL+27 : effective resolved shear stress

CFIX    13*NSLPTL+28                : volume fraction of void

CFIX    13*NSLPTL+29                : microscopic equivalent plastic strain

CFIX    13*NSLPTL+30 - NSTATV-9  : additional parameters users may need 
C                                  to characterize the constitutive law 
C                                  of a single crystal (if there are 
C                                  any).
C----   NSTATV-8               :  Plastic volume increase
C       NSTATV-6               :  number of iteration
C       NSTATV-5               :  cumulative equivalent plastic strain
C       NSTATV-4               :  determinant of deformation gradient
C       NSTATV-3               :  number of slip systems in the 1st set
C       NSTATV-2               :  number of slip systems in the 2nd set
C       NSTATV-1               :  number of slip systems in the 3rd set
C       NSTATV                 :  total number of slip systems in all sets

CFIX    NSTATV>=11*NSLPTL+35
C
C
C-----  Material constants PROPS:
C
C       PROPS(1) - PROPS(21) -- elastic constants for a general elastic
C                               anisotropic material
C
C            isotropic   : PROPS(i)=0  for  i>2
C                          PROPS(1) -- Young's modulus
C                          PROPS(2) -- Poisson's ratio
C
C            cubic       : PROPS(i)=0  for i>3
C                          PROPS(1) -- c11
C                          PROPS(2) -- c12
C                          PROPS(3) -- c44
C
C            orthotropic : PROPS(i)=0  for  i>9
C                          PROPS(1) - PROPS(9) are D1111, D1122, D2222,
C                          D1133, D2233, D3333, D1212, D1313, D2323, 
C                          respectively, which has the same definition 
C                          as ABAQUS for orthotropic materials
C                          (see *ELASTIC card)
C
C            anisotropic : PROPS(1) - PROPS(21) are D1111, D1122, 
C                          D2222, D1133, D2233, D3333, D1112, D2212, 
C                          D3312, D1212, D1113, D2213, D3313, D1213, 
C                          D1313, D1123, D2223, D3323, D1223, D1323, 
C                          D2323, respectively, which has the same 
C                          definition as ABAQUS for anisotropic 
C                          materials (see *ELASTIC card)
C
C
C       PROPS(25) - PROPS(56) -- parameters characterizing all slip 
C                                systems to be activated in a cubic 
C                                crystal
C
C            PROPS(25) -- number of sets of slip systems (maximum 3), 
C                         e.g. (110)[1-11] and (101)[11-1] are in the 
C                         same set of slip systems, (110)[1-11] and 
C                         (121)[1-11] belong to different sets of slip 
C                         systems
C                         (It must be a real number, e.g. 3., not 3 !)
C
C            PROPS(33) - PROPS(35) -- normal to a typical slip plane in
C                                     the first set of slip systems, 
C                                     e.g. (1 1 0)
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 0., not 1 1 0 !)
C            PROPS(36) - PROPS(38) -- a typical slip direction in the 
C                                     first set of slip systems, e.g. 
C                                     [1 1 1]
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 1., not 1 1 1 !)
C
C            PROPS(41) - PROPS(43) -- normal to a typical slip plane in
C                                     the second set of slip systems
C                                     (real numbers)
C            PROPS(44) - PROPS(46) -- a typical slip direction in the 
C                                     second set of slip systems
C                                     (real numbers)
C
C            PROPS(49) - PROPS(51) -- normal to a typical slip plane in
C                                     the third set of slip systems
C                                     (real numbers)
C            PROPS(52) - PROPS(54) -- a typical slip direction in the 
C                                     third set of slip systems
C                                     (real numbers)
C
C
C       PROPS(57) - PROPS(72) -- parameters characterizing the initial 
C                                orientation of a single crystal in 
C                                global system
C            The directions in global system and directions in local 
C            cubic crystal system of two nonparallel vectors are needed
C            to determine the crystal orientation.
C
C            PROPS(57) - PROPS(59) -- [p1 p2 p3], direction of first 
C                                     vector in local cubic crystal 
C                                     system, e.g. [1 1 0]
C                                     (They must be real numbers, e.g. 
C                                      1. 1. 0., not 1 1 0 !)
C            PROPS(60) - PROPS(62) -- [P1 P2 P3], direction of first 
C                                     vector in global system, e.g. 
C                                     [2. 1. 0.]
C                                     (It does not have to be a unit 
C                                      vector)
C
C            PROPS(65) - PROPS(67) -- direction of second vector in 
C                                     local cubic crystal system (real 
C                                     numbers)
C            PROPS(68) - PROPS(70) -- direction of second vector in 
C                                     global system
C
C
C       PROPS(73) - PROPS(96) -- parameters characterizing the visco-
C                                plastic constitutive law (shear 
C                                strain-rate vs. resolved shear 
C                                stress), e.g. a power-law relation
C
C            PROPS(73) - PROPS(80) -- parameters for the first set of 
C                                     slip systems
C            PROPS(81) - PROPS(88) -- parameters for the second set of 
C                                     slip systems
C            PROPS(89) - PROPS(96) -- parameters for the third set of 
C                                     slip systems
C
C
C       PROPS(97) - PROPS(144)-- parameters characterizing the self-
C                                and latent-hardening laws of slip 
C                                systems
C
C            PROPS(97) - PROPS(104)-- self-hardening parameters for the
C                                     first set of slip systems
C            PROPS(105)- PROPS(112)-- latent-hardening parameters for 
C                                     the first set of slip systems and
C                                     interaction with other sets of 
C                                     slip systems
C
C            PROPS(113)- PROPS(120)-- self-hardening parameters for the
C                                     second set of slip systems
C            PROPS(121)- PROPS(128)-- latent-hardening parameters for 
C                                     the second set of slip systems 
C                                     and interaction with other sets 
C                                     of slip systems
C
C            PROPS(129)- PROPS(136)-- self-hardening parameters for the
C                                     third set of slip systems
C            PROPS(137)- PROPS(144)-- latent-hardening parameters for 
C                                     the third set of slip systems and
C                                     interaction with other sets of
C                                     slip systems
C
C
C       PROPS(145)- PROPS(152)-- parameters characterizing forward time
C                                integration scheme and finite 
C                                deformation
C
C            PROPS(145) -- parameter theta controlling the implicit 
C                          integration, which is between 0 and 1
C                          0.  : explicit integration
C                          0.5 : recommended value
C                          1.  : fully implicit integration
C
C            PROPS(146) -- parameter NLGEOM controlling whether the 
C                          effect of finite rotation and finite strain 
C                          of crystal is considered,
C                          0.        : small deformation theory
C                          otherwise : theory of finite rotation and 
C                                      finite strain
C
C
C       PROPS(153)- PROPS(160)-- parameters characterizing iteration 
C                                method
C
C            PROPS(153) -- parameter ITRATN controlling whether the 
C                          iteration method is used, 
C                          0.        : no iteration
C                          otherwise : iteration
C
C            PROPS(154) -- maximum number of iteration ITRMAX 
C
C            PROPS(156) -- absolute error of ΔραSSD in slip 
C                          systems DRHOERR
C            PROPS(157) -- maximum number of iteration for calculating ΔραSSD
C

      PI=ACOS(-1.0)

C-----  Elastic matrix in local cubic crystal system: DLOCAL
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.
         END DO
      END DO

      C11=PROPS(1)+PROPS(4)*PROPS(10)
      C12=PROPS(2)+PROPS(5)*PROPS(10)
      C44=PROPS(3)+PROPS(6)*PROPS(10)

      DO I=1,3
         DO J=1,3
            IF (I.EQ.J) THEN
               DLOCAL(I,J)=C11
               DLOCAL(I+3,J+3)=C44
            ELSE
               DLOCAL(I,J)=C12
            END IF
         END DO
      END DO

      SHEARMD=(1./5.*(C11-C12+3.*C44)
     2         +5.*C44*(C11-C12)/(4.*C44+3.*(C11-C12)))/2.



      CALL ROTATION (PROPS(57), ROTATE)



C-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL 
C     to global elastic matrix D
C
      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

C-----  Elastic matrix in global system: D
C     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose
C
      DO J=1,6
         DO I=1,6
            D(I,J)=0.
         END DO
      END DO

      DO J=1,6
         DO I=1,J

            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO

            D(J,I)=D(I,J)

         END DO
      END DO

C-----  Total number of sets of slip systems: NSET
      NSET=NINT(PROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         STOP
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*) 
     2     '***ERROR - more than three sets of slip systems'
         STOP
      END IF

C-----  Implicit integration parameter: THETA
      THETA=PROPS(145)


C-----  Finite deformation ?
C-----  NLGEOM = 0,   small deformation theory
C       otherwise, theory of finite rotation and finite strain, Users 
C     must declare "NLGEOM" in the input file, at the *STEP card
C
      IF (PROPS(146).EQ.0.) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF



C-----  Iteration?
C-----  ITRATN = 0, no iteration
C       otherwise, iteration (solving increments of stresses and 
C     solution dependent state variables)
C
      IF (PROPS(153).EQ.0.) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF

      ITRMAX=NINT(PROPS(154))
      GAMERR=PROPS(155)

      NITRTN=-1

      DO I=1,NTENS
         DSOLD(I)=0.
         DELATSOD(I)=0.
         DPLATSOD(I)=0.
      END DO

      DO J=1,ND
         DGAMOD(J)=0.
         DTAUOD(J)=0.
         DRHOOD(J)=0.
         DO I=1,3
            DSPNRO(I,J)=0.
            DSPDRO(I,J)=0.
         END DO
      END DO

      DEQPLATSOD=0.
      DMEQPLATSOD=0.
      DFVOIDOD=0.

C-----  Increment of spin associated with the material element: DSPIN
C     (only needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I)
               TRM0(I,J)=DROT(J,I)
            END DO

            TERM(J,J)=TERM(J,J)+1.D0
            TRM0(J,J)=TRM0(J,J)-1.D0
         END DO

         PRTSTR='DSPIN'
         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP,PRTSTR, IDPNEWDT)
         IF (IDPNEWDT.EQ.1) THEN
             PNEWDT=0.5
             GO TO 2002
         END IF

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         DSPIN(1)=TRM0(2,1)-TRM0(1,2)
         DSPIN(2)=TRM0(1,3)-TRM0(3,1)
         DSPIN(3)=TRM0(3,2)-TRM0(2,3)

      END IF

C-----  Increment of dilatational strain: DEV
      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO


C-----  Iteration starts (only when iteration method is used)
2001  CONTINUE

C-----  Parameter NITRTN: number of iterations
C       NITRTN = 0 --- no-iteration solution
C
      NITRTN=NITRTN+1

C-----  Check whether the current stress state is the initial state
      IF (STATEV(1).EQ.0.) THEN

C-----  Initial state
C
C-----  Generating the following parameters and variables at initial 
C     state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Unit vectors in initial slip directions SLPDIR
C          Unit normals to initial slip planes SLPNOR


         NSLPTL=0
         DO I=1,NSET
            ISPNOR(1)=NINT(PROPS(25+8*I))
            ISPNOR(2)=NINT(PROPS(26+8*I))
            ISPNOR(3)=NINT(PROPS(27+8*I))

            ISPDIR(1)=NINT(PROPS(28+8*I))
            ISPDIR(2)=NINT(PROPS(29+8*I))
            ISPDIR(3)=NINT(PROPS(30+8*I))


            CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), SLPDIR(1,NSLPTL+1), 
     2                    SLPNOR(1,NSLPTL+1), ROTATE)

            NSLPTL=NSLPTL+NSLIP(I)
         END DO

         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*) 
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            STOP
         END IF



C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
            
            SLPDEFNOR(1,J)=SLPNOR(1,J)*SLPNOR(1,J)
            SLPDEFNOR(2,J)=SLPNOR(2,J)*SLPNOR(2,J)
            SLPDEFNOR(3,J)=SLPNOR(3,J)*SLPNOR(3,J)
            SLPDEFNOR(4,J)=SLPNOR(1,J)*SLPNOR(2,J)
     2                    +SLPNOR(2,J)*SLPNOR(1,J)
            SLPDEFNOR(5,J)=SLPNOR(1,J)*SLPNOR(3,J)
     2                    +SLPNOR(3,J)*SLPNOR(1,J)
            SLPDEFNOR(6,J)=SLPNOR(2,J)*SLPNOR(3,J)
     2                    +SLPNOR(3,J)*SLPNOR(2,J)
         END DO



C-----  Initial value of state variables: unit normal to a slip plane 
C       and unit vector in a slip direction

         STATEV(NSTATV)=FLOAT(NSLPTL)
         DO I=1,NSET
            STATEV(NSTATV-4+I)=FLOAT(NSLIP(I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=SLPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=SLPDIR(I,J)
            END DO
         END DO


C-----  Initial value of the resolved shear stress in slip systems
         DO I=1,NSLPTL
            TERM1=0.

            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO

            STATEV(2*NSLPTL+I)=TERM1
         END DO

C-----  Initial value of the effective resolved shear stress in slip systems
         DO I=1,NSLPTL
            STATEV(12*NSLPTL+27+I)=0.
         END DO


C-----  Initial value of the effective resolved normal stress in slip systems
         DO I=1,NSLPTL
            STATEV(10*NSLPTL+I)=PROPS(10)
         END DO


C-----  Initial value of shear strain in slip systems
CFIX--  Initial value of cumulative shear strain in each slip systems
        DO I=1,NSLPTL
           STATEV(NSLPTL+I)=0.
           STATEV(9*NSLPTL+I)=0.
        END DO


C-----  Initial value of elastic, plastic, total strain components
        DO I=1,NTENS
           STATEV(11*NSLPTL+I)=0.
           STATEV(11*NSLPTL+6+I)=0.
           STATEV(11*NSLPTL+12+I)=0.
        END DO


C-----  Initial value of dislocation density in each slip system
       ID=0
       DO I=1,NSET
          DO J=1,NSLIP(I)
             ID=ID+1
             STATEV(11*NSLPTL+27+ID)=PROPS(96+(I-1)*9)
          END DO
       END DO


C----  Coefficient matrix: SLIPDM
      CALL SLIPSYSCOE (SLIPDM,ND,NSLPTL,PROPS(73))


C-----  Initial value of the current strength for all slip systems
       ID=0
       DO I=1,NSET
          DO J=1,NSLIP(I)
             ID=ID+1
             VAR1=0.
             DO K=1,NSLPTL
                VAR1=VAR1+SLIPDM(ID,K)*STATEV(11*NSLPTL+27+K)
             END DO
             STATEV(ID)=(PROPS(93+(I-1)*9)
     2                   +PROPS(92)*SHEARMD
     3                    *PROPS(95+(I-1)*9)*(VAR1**0.5))
     4                  *(1.+PROPS(9)*STATEV(10*NSLPTL+ID))
          END DO
       END DO

C-----  Initial Weibull distribution
        ! open(unit=99,file='D:\FuJiaqi\W4-Swelling\Weibull\k3.dat')
        ! DO I = 1,50000
            ! READ(99,*) Weibull(I)
        ! END DO
        ! close(unit=99)

C-----  Initial value of equivalent plastic strain
        STATEV(11*NSLPTL+19)=0.
        
C-----  Initial value of cumulative equivalent plastic strain
        STATEV(NSTATV-5)=0.

C-----  Initial plastic volume increase
        STATEV(NSTATV-8)=0.

C-----  Initial equivalent strain
        STATEV(11*NSLPTL+27)=0.

C-----  Initial volume fraction of void
        ! RANDUM_NUM = int(50000*RANDUM_NUM)
        ! IF(RANDUM_NUM.eq.0)THEN
            ! RANDUM_NUM = 50000
        ! END IF
        if (cmname(1:2).eq.'GB') then
          STATEV(13*NSLPTL+28)=PROPS(20)*2.
        else
          STATEV(13*NSLPTL+28)=PROPS(20)*0.903994
        end if
        IF(STATEV(13*NSLPTL+28).lt.0)THEN
        STATEV(13*NSLPTL+28) = 0.
        END IF
C-----  Initial value of microscopic equivalent plastic strain
        STATEV(13*NSLPTL+29)=0.

CFIX    initial total dislocation density on all slip systems    
         STATEV(11*NSLPTL+20)=0.
         DO I=1,NSLPTL
            STATEV(11*NSLPTL+20)=STATEV(11*NSLPTL+20)
     2                           +STATEV(11*NSLPTL+27+I)
         END DO

C-----  Initial value of dislocation density in slip system
       DO I=1,NSET
          STATEV(NSTATV-8-I)=PROPS(96+(I-1)*9)
       END DO

C-----  Initial grain size
       STATEV(NSTATV-8-NSET-1)=PROPS(137)


      ELSE


C-----  Current stress state
C
C-----  Copying from the array of state variables STATEV the following
C          parameters and variables at current stress state:
C          Total number of slip systems in all the sets NSLPTL
C          Number of slip systems in each set NSLIP
C          Current slip directions SLPDIR
C          Normals to current slip planes SLPNOR
C
         NSLPTL=NINT(STATEV(NSTATV))
         DO I=1,NSET
            NSLIP(I)=NINT(STATEV(NSTATV-4+I))
         END DO


         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               SLPNOR(I,J)=STATEV(IDNOR)

               IDDIR=IDDIR+1
               SLPDIR(I,J)=STATEV(IDDIR)
            END DO
         END DO



C-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
            
            SLPDEFNOR(1,J)=SLPNOR(1,J)*SLPNOR(1,J)
            SLPDEFNOR(2,J)=SLPNOR(2,J)*SLPNOR(2,J)
            SLPDEFNOR(3,J)=SLPNOR(3,J)*SLPNOR(3,J)
            SLPDEFNOR(4,J)=SLPNOR(1,J)*SLPNOR(2,J)
     2                    +SLPNOR(2,J)*SLPNOR(1,J)
            SLPDEFNOR(5,J)=SLPNOR(1,J)*SLPNOR(3,J)
     2                    +SLPNOR(3,J)*SLPNOR(1,J)
            SLPDEFNOR(6,J)=SLPNOR(2,J)*SLPNOR(3,J)
     2                    +SLPNOR(3,J)*SLPNOR(2,J)
            
         END DO

      END IF



C-----  Slip spin tensor: SLPSPN (only needed for finite rotation)
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF



C-----  Double dot product of elastic moduli tensor with the slip 
C     deformation tensor (Schmid factors) plus, only for finite 
C     rotation, the dot product of slip spin tensor with the stress: 
C     DDEMSD (Pα:Ce)
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF


      ID=0
      DO P=1,NSET
         DO Q=1,NSLIP(P)

            ID=ID+1
            DO I=1,6
               DDEMSDNOR(I,ID)=0.
               DO K=1,6
                  DDEMSDNOR(I,ID)=DDEMSDNOR(I,ID)
     2                            +D(K,I)*SLPDEFNOR(K,ID)
               END DO
            END DO

         END DO
      END DO

C     calculate Mises stress, hydrostatic pressure, deviatoric stress
      HYDROSTRESS=(STRESS(1)+STRESS(2)+STRESS(3))/3.
      !!! 静水压力
      DEVSIG(1,1)=STRESS(1)-HYDROSTRESS
      DEVSIG(2,2)=STRESS(2)-HYDROSTRESS
      DEVSIG(3,3)=STRESS(3)-HYDROSTRESS
      DEVSIG(1,2)=STRESS(4)
      DEVSIG(1,3)=STRESS(5)
      DEVSIG(2,3)=STRESS(6)
      DEVSIG(2,1)=STRESS(4)
      DEVSIG(3,1)=STRESS(5)
      DEVSIG(3,2)=STRESS(6)
      !!! 偏应力的分量
      EQSTRESS=0.
      DO I=1,3
         DO J=1,3
            EQSTRESS=EQSTRESS+DEVSIG(I,J)*DEVSIG(I,J)
         END DO
      END DO

      EQSTRESS=(3./2.*EQSTRESS)**0.5
      !!! von-Mises应力

C---- calculate effective void fraction
      IF (STATEV(13*(NSLPTL)+28).GE.(1./PROPS(12))) THEN
          !!! STATEV(13*(NSLPTL)+28) 孔隙率
          !!! PROPS(12) q1
          FVOIDEFF=1./PROPS(12)
      ELSE IF (STATEV(13*(NSLPTL)+28).LE.0.) THEN
          FVOIDEFF=0.
      ELSE
          FVOIDEFF=STATEV(13*(NSLPTL)+28)
      END IF
      !!! 这一段是对孔隙率的标准化

      DO I=1,NSLPTL
         IDTAUEFF(I)=0
      END DO
      !!! 是否存在等效分切应力

      IF (EQSTRESS.EQ.0.) THEN
         DO I=1,NSLPTL
            TAUSPEFF(I)=0.
            !!! 如果等效塑性应力为0，则等效分切应力为0
         END DO
      ELSE
         DO I=1,NSLPTL
            IF (STATEV(2*NSLPTL+I).EQ.0.) THEN
                TAUSPEFF(I)=0.
                !!! 如果分切应力是0
                !!! 则等效分切应力为0
            ELSE
                ID=-1
9001            CONTINUE
                ID=ID+1
                IF (ID.EQ.0) THEN
                   TAUSPEFF(I)=abs(STATEV(2*NSLPTL+I))
                   !!! 等效分切应力的初始值
                   !!! 就是分切应力
                ELSE
                   !!! 牛顿迭代法
                   !!! PROPS(11) lambda
                   !!! PROPS(12) q1
                   !!! PROPS(13) q2
                   VAR1=(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     2                 +2.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
     3                    *DSIGN(1.0,STATEV(2*NSLPTL+I))
     4                    /(TAUSPEFF(I))**2.
     2                 +2./45.*PROPS(11)*FVOIDEFF
     3                        *(EQSTRESS/TAUSPEFF(I))**2.
     4                 +2.*PROPS(12)*FVOIDEFF
     5                   *(1.+0.5*PROPS(13)**2.*(
     6                   3./10.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**2.+
     7                   3./20.*(HYDROSTRESS/TAUSPEFF(I))**2.))
     7                 -1.-(PROPS(12)*FVOIDEFF)**2.

                   VAR3=4.*1./TAUSPEFF(I)
                   VAR3=VAR3*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
                   VAR3=VAR3*DSIGN(1.0,STATEV(2*NSLPTL+I))
                   VAR3=VAR3/TAUSPEFF(I)**2.

                   VAR2=-2.*1./TAUSPEFF(I)
     2                     *(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     3                  -VAR3
     3                  -4./45.*PROPS(11)*FVOIDEFF
     4                    *1./TAUSPEFF(I)*(EQSTRESS/TAUSPEFF(I))**2.
     5                  -2.*PROPS(12)*FVOIDEFF*0.5*PROPS(13)**2.*
     6                   (3./5.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**3.
     7                   +3./10.*HYDROSTRESS**2./TAUSPEFF(I)**3.)

                   TAUSPEFF(I)=TAUSPEFF(I)-VAR1/VAR2
                END IF
                !!! RESIDU 看起来是残差
                !!! 实际上就是屈服势函数和0之间的差值
                RESIDU=(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     2                 +2.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
     3                    *DSIGN(1.0,STATEV(2*NSLPTL+I))
     4                    /(TAUSPEFF(I))**2.
     2                 +2./45.*PROPS(11)*FVOIDEFF
     3                        *(EQSTRESS/TAUSPEFF(I))**2.
     4                 +2.*PROPS(12)*FVOIDEFF
     5                   *(1.+0.5*PROPS(13)**2.*(
     6                   3./10.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**2.+
     7                   3./20.*(HYDROSTRESS/TAUSPEFF(I))**2.))
     7                 -1.-(PROPS(12)*FVOIDEFF)**2.

                IF (ISNAN(RESIDU).OR.(ABS(RESIDU).GT.(10.**10.))) THEN
                   TAUSPEFF(I)=STATEV(2*NSLPTL+I)
                   IDTAUEFF(I)=1
                   !!! 如果说残差不收敛
                   !!! 则等效分切应力就是分切应力
                   GO TO 9002
                ELSE IF (ABS(RESIDU).GT.(10.**-6.)) THEN
                   !!! 如果残差较大
                   !!! 就继续迭代
                   !!! 看起来ID就是迭代的次数
                   GO TO 9001
                ELSE IF (ABS(RESIDU).LE.(10.**-6.)) THEN
                   !!! 如果残差比1e-6小了
                   !!! 那么就计算得到了正确的等效分切应力
                   GO TO 9002
                   IDTAUEFF(I)=2
                   !!! 如果说表示符为2，则说明得到了正确的等效分切应力
                END IF

            END IF

            TAUSPEFF(I) = TAUSPEFF(I)*DSIGN(1.0,STATEV(2*NSLPTL+I))

9002        CONTINUE
         END DO

      END IF

      SLPDEFEFF33=0.

      IF (EQSTRESS.EQ.0.) THEN
          DO I=1,NSLPTL
             DO J=1,3
                DO K=1,3
                   SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                END DO
             END DO
          END DO
      ELSE
          DO I=1,NSLPTL
             IF (TAUSPEFF(I).EQ.0.) THEN
                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                       !!! 如果等效分切应力等于0
                       !!! 滑移系张量为施密特张量
                    END DO
                 END DO
             ELSE IF (IDTAUEFF(I).EQ.1) THEN
                 DO J=1,3
                    DO K=1,3
                       !!! 如果不能得到等效分切应力
                       !!! 滑移系张量为施密特张量
                       SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                    END DO
                 END DO
             ELSE
                 VAR1=-2.*1./TAUSPEFF(I)
     2                     *(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     3                  -4.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)*
     4                     *DSIGN(1.0,STATEV(2*NSLPTL+I))
     5                     /TAUSPEFF(I)**3.
     3                  -4./45.*PROPS(11)*FVOIDEFF
     4                    *1./TAUSPEFF(I)*(EQSTRESS/TAUSPEFF(I))**2.
     5                  -2.*PROPS(12)*FVOIDEFF*0.5*PROPS(12)**2.*
     6                   (3./5.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**3.
     7                   +3./10.*HYDROSTRESS**2./TAUSPEFF(I)**3.)

                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=2.*STATEV(2*NSLPTL+I)
     2                                      /TAUSPEFF(I)**2.
     3                                     *SLPDIR(J,I)*SLPNOR(K,I)
     4                                   +2.*FVOIDEFF*PROPS(17)
     5                                      *DSIGN(1.0,STATEV(2*NSLPTL+I))
     6                                      /TAUSPEFF(I)**2.
     7                                     *SLPDIR(J,I)*SLPNOR(K,I)
     4                                   +2./15.*PROPS(11)*FVOIDEFF
     5                                  /TAUSPEFF(I)**2.*DEVSIG(J,K)
                    END DO
                    SLPDEFEFF33(J,J,I)=SLPDEFEFF33(J,J,I)
     5                   +2.*PROPS(12)*FVOIDEFF*0.5*PROPS(13)**2.
     6                    /10.*(PROPS(17)+HYDROSTRESS)
     7                    /TAUSPEFF(I)**2.
                 END DO
     
                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=-SLPDEFEFF33(J,K,I)/VAR1
                       !!! 这个就是考虑孔隙率效应的滑移系张量
                    END DO
                 END DO

             END IF
          END DO
      END IF


      SLPDEFNOREFF33=0.
      !!! 将滑移系张量正则化
      DO I=1,NSLPTL
         DO J=1,3
            VAR1=SLPDEFEFF33(J,1,I)**2.+SLPDEFEFF33(J,2,I)**2.
     2           +SLPDEFEFF33(J,3,I)**2.
            IF (VAR1.NE.0.) THEN
               DO K=1,3
                  DO L=1,3
                     SLPDEFNOREFF33(K,L,I)=
     2               SLPDEFEFF33(J,K,I)*SLPDEFEFF33(J,L,I)/VAR1
                  END DO
               END DO
               GO TO 7001
            END IF
         END DO
7001     CONTINUE
      END DO

      !!! 以下为滑移系张量的对称部分

      SLPDEFEFF=0.
      DO I=1,NSLPTL
         SLPDEFEFF(1,I)=SLPDEFEFF33(1,1,I)
         SLPDEFEFF(2,I)=SLPDEFEFF33(2,2,I)
         SLPDEFEFF(3,I)=SLPDEFEFF33(3,3,I)
         SLPDEFEFF(4,I)=SLPDEFEFF33(1,2,I)+SLPDEFEFF33(2,1,I)
         SLPDEFEFF(5,I)=SLPDEFEFF33(1,3,I)+SLPDEFEFF33(3,1,I)
         SLPDEFEFF(6,I)=SLPDEFEFF33(2,3,I)+SLPDEFEFF33(3,2,I)
      END DO
      

      SLPDEFNOREFF=0.
      DO I=1,NSLPTL
         SLPDEFNOREFF(1,I)=SLPDEFNOREFF33(1,1,I)
         SLPDEFNOREFF(2,I)=SLPDEFNOREFF33(2,2,I)
         SLPDEFNOREFF(3,I)=SLPDEFNOREFF33(3,3,I)
         SLPDEFNOREFF(4,I)=SLPDEFNOREFF33(1,2,I)
     2                    +SLPDEFNOREFF33(2,1,I)
         SLPDEFNOREFF(5,I)=SLPDEFNOREFF33(1,3,I)
     2                    +SLPDEFNOREFF33(3,1,I)
         SLPDEFNOREFF(6,I)=SLPDEFNOREFF33(2,3,I)
     2                    +SLPDEFNOREFF33(3,2,I)
      END DO

      !!! 以下为滑移系张量的非对称部分（旋转张量）

      IF (NLGEOM.NE.0) THEN
          DO I=1,NSLPTL
             SLPSPNEFF(1,I)=0.5*(SLPDEFEFF33(1,2,I)-
     2                           SLPDEFEFF33(2,1,I))
             SLPSPNEFF(2,I)=0.5*(SLPDEFEFF33(3,1,I)-
     2                           SLPDEFEFF33(1,3,I))
             SLPSPNEFF(3,I)=0.5*(SLPDEFEFF33(2,3,I)-
     3                           SLPDEFEFF33(3,2,I))
          END DO
      END IF




C-----  Double dot product of elastic moduli tensor with the effective slip 
C     deformation tensor (effective Schmid factors) plus, only for finite 
C     rotation, the dot product of slip spin tensor with the stress: 
C     DDEMSDEFF (Pαeff:Ce)
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSDEFF(I,J)=0.
            DO K=1,6
               DDEMSDEFF(I,J)=DDEMSDEFF(I,J)+D(K,I)*SLPDEFEFF(K,J)
            END DO
         END DO
      END DO

      DO J=1,NSLPTL
         DO I=1,6
            DDEMSDEFF0(I,J)=0.
            DO K=1,6
               DDEMSDEFF0(I,J)=DDEMSDEFF0(I,J)+D(K,I)*SLPDEFEFF(K,J)
            END DO
         END DO
      END DO



      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSDEFF(4,J)=DDEMSDEFF(4,J)-SLPSPNEFF(1,J)*STRESS(1)
            DDEMSDEFF(5,J)=DDEMSDEFF(5,J)+SLPSPNEFF(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSDEFF(4,J)=DDEMSDEFF(4,J)+SLPSPNEFF(1,J)*STRESS(2)
               DDEMSDEFF(6,J)=DDEMSDEFF(6,J)-SLPSPNEFF(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSDEFF(5,J)=DDEMSDEFF(5,J)-SLPSPNEFF(2,J)*STRESS(3)
               DDEMSDEFF(6,J)=DDEMSDEFF(6,J)+SLPSPNEFF(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
              DDEMSDEFF(1,J)=DDEMSDEFF(1,J)+SLPSPNEFF(1,J)*STRESS(NDI+1)
              DDEMSDEFF(2,J)=DDEMSDEFF(2,J)-SLPSPNEFF(1,J)*STRESS(NDI+1)
              DDEMSDEFF(5,J)=DDEMSDEFF(5,J)-SLPSPNEFF(3,J)*STRESS(NDI+1)
              DDEMSDEFF(6,J)=DDEMSDEFF(6,J)+SLPSPNEFF(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
              DDEMSDEFF(1,J)=DDEMSDEFF(1,J)-SLPSPNEFF(2,J)*STRESS(NDI+2)
              DDEMSDEFF(3,J)=DDEMSDEFF(3,J)+SLPSPNEFF(2,J)*STRESS(NDI+2)
              DDEMSDEFF(4,J)=DDEMSDEFF(4,J)+SLPSPNEFF(3,J)*STRESS(NDI+2)
              DDEMSDEFF(6,J)=DDEMSDEFF(6,J)-SLPSPNEFF(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
              DDEMSDEFF(2,J)=DDEMSDEFF(2,J)+SLPSPNEFF(3,J)*STRESS(NDI+3)
              DDEMSDEFF(3,J)=DDEMSDEFF(3,J)-SLPSPNEFF(3,J)*STRESS(NDI+3)
              DDEMSDEFF(4,J)=DDEMSDEFF(4,J)-SLPSPNEFF(2,J)*STRESS(NDI+3)
              DDEMSDEFF(5,J)=DDEMSDEFF(5,J)+SLPSPNEFF(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF

      DO J=1,NSLPTL
         DO I=1,6
            DDEMSDNOREFF(I,J)=0.
            DO K=1,6
               DDEMSDNOREFF(I,J)=DDEMSDNOREFF(I,J)
     2                           +D(K,I)*SLPDEFNOREFF(K,J)
            END DO
         END DO
      END DO


C----  Coefficient matrix: SLIPDM
      CALL SLIPSYSCOE (SLIPDM,ND,NSLPTL,PROPS(73))

C----  Coefficient matrix: SLIPGM
      CALL SLIPSYSCOE (SLIPGM,ND,NSLPTL,PROPS(79))

  
C-----  Shear strain-rate in a slip system at the start of increment 
C     FSLIP, and its derivative DFDXTAU, DFDXRHO, DFDXSIG
       
       IF (EQSTRESS.EQ.0.) THEN
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)

                ID=ID+1
                FSLIP(ID)=0.
                DFDXTAUEFF(ID)=0.
                DTAUEFFDTAU(ID)=0.
                DTAUEFFDFVOID(ID)=0.
                DFDXSIGEFF(ID)=0.

                DO K=1,NSLPTL
                   DFDXRHO(ID,K)=0.
                END DO

             END DO
          END DO

       ELSE
          !!! 有效孔洞体积分数和孔洞体积分数之间的比例
          IF (STATEV(13*(NSLPTL)+28).GE.(1./PROPS(12))) THEN
              DFVOIDEFFDFVOID=0.
          ELSE IF (STATEV(13*(NSLPTL)+28).LE.0.) THEN
              DFVOIDEFFDFVOID=0.
          ELSE
              DFVOIDEFFDFVOID=1.
          END IF

          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)

                ID=ID+1
                IF (TAUSPEFF(ID).EQ.0.) THEN
                    
                    FSLIP(ID)=0.
                    DFDXTAUEFF(ID)=0.
                    DTAUEFFDTAU(ID)=0.
                    DTAUEFFDFVOID(ID)=0.
                    DFDXSIGEFF(ID)=0.
                    
                    DO K=1,NSLPTL
                       DFDXRHO(ID,K)=0.
                    END DO
                    
                ELSE IF (IDTAUEFF(ID).EQ.1) THEN

C                    FSLIP(ID)=0.
C                    DFDXTAUEFF(ID)=0.
C                    DTAUEFFDTAU(ID)=0.
C                    DTAUEFFDFVOID(ID)=0.
C                    DFDXSIGEFF(ID)=0.
C                    
C                    DO K=1,NSLPTL
C                       DFDXRHO(ID,K)=0.
C                    END DO

                    VAR1=0.
                    DO K=1,NSLPTL
                       VAR1=VAR1+SLIPDM(ID,K)*STATEV(11*NSLPTL+27+K)
                    END DO
                    !!! STATEV(11*NSLPTL+27+K) 各个滑移系上的位错密度
                    !!! SLIPDM(ID,K) 位错的相互作用强度

                    VAR2=(PROPS(93+(I-1)*9)
     2                    +PROPS(92)*SHEARMD
     3                     *PROPS(95+(I-1)*9)*(VAR1**0.5))
     4                   *(1.+PROPS(9)*STATEV(10*NSLPTL+ID))
                    !!! (PROPS(93+(I-1)*9) 常数的滑移阻力
                    !!! PROPS(92) 位错强度c
                    !!! PROPS(95+(I-1)*9) 柏氏矢量
                    !!! PROPS(9)这个9是师兄为了考虑静水压力作用
                    !!! 引入的系数，一定要注意，时刻把它设为0

                    VAR3=PROPS(93+(I-1)*9)
     2                   +PROPS(92)*SHEARMD
     3                    *PROPS(95+(I-1)*9)*(VAR1**0.5)
                    !!! 去除静水压力作用的VAR2
                    FSLIP(ID)=PROPS(99+(I-1)*9)
     2                        *(ABS(TAUSPEFF(ID)/VAR2))
     3                          **PROPS(100+(I-1)*9)
     4                        *DSIGN(1.D0,TAUSPEFF(ID))
                    !!! PROPS(99+(I-1)*9) 参考应变率
                    !!! TAUSPEFF(ID) 等效分切应力
                    !!! PROPS(100+(I-1)*9) 应变率敏感性因子
                    DFDXTAUEFF(ID)=PROPS(100+(I-1)*9)*PROPS(99+(I-1)*9)
     2                   /VAR2*(ABS(TAUSPEFF(ID)/VAR2))
     3                   **(PROPS(100+(I-1)*9)-1.)
                    !!! DFDXTAUEFF(ID) 滑移率对等效分切应力的导数
                    DTAUEFFDTAU(ID)=1.
                    DTAUEFFDFVOID(ID)=0.

                    DFDXSIGEFF(ID)=-PROPS(100+(I-1)*9)*PROPS(9)
     2                    /(1.+PROPS(9)*STATEV(10*NSLPTL+ID))*FSLIP(ID)

                    !!! 滑移率对位错密度的导数
                    VAR2=-PROPS(100+(I-1)*9)*PROPS(92)*SHEARMD
     2                    *PROPS(95+(I-1)*9)
     3                    /(2.*VAR3*VAR1**0.5)*FSLIP(ID)
                    !!! 滑移率对位错密度的导数(乘上了交叉矩阵)
                    DO K=1,NSLPTL
                       DFDXRHO(ID,K)=SLIPDM(ID,K)*VAR2
                    END DO
                    
                ELSE
                
                    VAR1=0.
                    DO K=1,NSLPTL
                       VAR1=VAR1+SLIPDM(ID,K)*STATEV(11*NSLPTL+27+K)
                    END DO
                    !!! STATEV(11*NSLPTL+27+K) 各个滑移系上的位错密度
                    !!! SLIPDM(ID,K) 位错的相互作用强度
                    
                    VAR2=(PROPS(93+(I-1)*9)
     2                    +PROPS(92)*SHEARMD
     3                    *PROPS(95+(I-1)*9)*(VAR1**0.5))
     4                    *(1.+PROPS(9)*STATEV(10*NSLPTL+ID))
                    
                    VAR3=PROPS(93+(I-1)*9)
     2                   +PROPS(92)*SHEARMD
     3                    *PROPS(95+(I-1)*9)*(VAR1**0.5)
                    
                    FSLIP(ID)=PROPS(99+(I-1)*9)
     2                        *(ABS(TAUSPEFF(ID)/VAR2))
     3                          **PROPS(100+(I-1)*9)
     4                        *DSIGN(1.D0,TAUSPEFF(ID))
         
                    DFDXTAUEFF(ID)=PROPS(100+(I-1)*9)*PROPS(99+(I-1)*9)
     2                          /VAR2*(ABS(TAUSPEFF(ID)/VAR2))
     3                          **(PROPS(100+(I-1)*9)-1.)
         
                    !!! phi对分切应力的导数
                    DPHIDTAU(ID)=2.*STATEV(2*NSLPTL+ID)
     2                           /(TAUSPEFF(ID)**2.)
                    !!! phi对等效分切应力的导数
                    DPHIDTAUEFF(ID)=-2.*1./TAUSPEFF(ID)
     2                     *(STATEV(2*NSLPTL+ID)/TAUSPEFF(ID))**2.
     3                  -4.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+ID)*
     4                     *DSIGN(1.0,STATEV(2*NSLPTL+ID))
     5                     /TAUSPEFF(ID)**3.
     3                  -4./45.*PROPS(11)*FVOIDEFF
     4                    *1./TAUSPEFF(ID)*(EQSTRESS/TAUSPEFF(ID))**2.
     5                  -2.*PROPS(12)*FVOIDEFF*0.5*PROPS(12)**2.*
     6                   (3./5.*PROPS(17)*HYDROSTRESS/TAUSPEFF(ID)**3.
     7                   +3./10.*HYDROSTRESS**2./TAUSPEFF(ID)**3.)

                    !!! phi对等效孔隙率的导数
                    DPHIDFVOIDEFF(ID)=2./45.*PROPS(11)
     2                                *(EQSTRESS/TAUSPEFF(ID))**2.
     3                   +2.*PROPS(17)*STATEV(2*NSLPTL+ID)
     4                      *DSIGN(1.0,STATEV(2*NSLPTL+ID))
     5                      /(TAUSPEFF(ID))**2.
     3                   +2.*PROPS(12)*(1.+0.5*PROPS(13)**2.*(
     4                    3./10.*PROPS(17)*HYDROSTRESS/TAUSPEFF(ID)**2.+
     5                    3./20.*(HYDROSTRESS/TAUSPEFF(ID))**2.))
     5                   -2.*PROPS(12)**2.*FVOIDEFF
         
                    !!! 等效分切应力对分切应力的导数
                    DTAUEFFDTAU(ID)=-DPHIDTAU(ID)/DPHIDTAUEFF(ID)
                    !!! 等效分切应力对孔隙率的导数
                    DTAUEFFDFVOID(ID)=-DPHIDFVOIDEFF(ID)*DFVOIDEFFDFVOID
     2                                 /DPHIDTAUEFF(ID)
                    !!! 由于不考虑压力(确保PROPS(9)=0)
                    DFDXSIGEFF(ID)=-PROPS(100+(I-1)*9)*PROPS(9)
     2                    /(1.+PROPS(9)*STATEV(10*NSLPTL+ID))*FSLIP(ID)
                    !!! 滑移率对位错密度的导数
                    VAR2=-PROPS(100+(I-1)*9)*PROPS(92)*SHEARMD
     2                    *PROPS(95+(I-1)*9)
     3                    /(2.*VAR3*VAR1**0.5)*FSLIP(ID)
                    
                    DO K=1,NSLPTL
                       DFDXRHO(ID,K)=SLIPDM(ID,K)*VAR2
                    END DO
                
                END IF


                
             END DO
          END DO

       END IF






C-----  LU decomposition to solve the increment of shear strain in a 
C     slip system
C
      IF (NITRTN.EQ.0) THEN

          TERM1=THETA*DTIME
          ID=0
          DO P=1,NSET
             DO Q=1,NSLIP(P)
                ID=ID+1
                TERM2=DFDXTAUEFF(ID)
                !!! 滑移率对等效分切应力的导数
                TERM3=DTAUEFFDTAU(ID)
                !!! 等效分切应力对分切应力的导数
                TERM4=DTAUEFFDFVOID(ID)
                !!! 等效分切应力对孔隙率的导数
                TERM5=DFDXSIGEFF(ID)
                !!! 等效分切应力对(Sigeff)的导数

                DO J=1,NSLPTL

                   TERM6=0.
                   DO K=1,6
                      TERM6=TERM6+DDEMSDEFF0(K,ID)
     2                *(1.-STATEV(13*NSLPTL+28))
     3                *(SLPDEFEFF(K,J)+PROPS(8)*SLPDEFNOREFF(K,J)
     4                                      *DSIGN(1.D0,FSLIP(J)))
                   !!! 这个是师兄Word的公式(42)
                   !!! PROPS(8) beta 需要设为0
                   END DO

                   TERM7=0.
                   DO K=1,6
                      TERM7=TERM7+
     2                DDEMSDNOREFF(K,ID)*(1.-STATEV(13*NSLPTL+28))
     3                *(SLPDEFEFF(K,J)+PROPS(8)*SLPDEFNOREFF(K,J)
     4                                 *DSIGN(1.D0,FSLIP(J)))
                   !!! 这个是师兄Word的公式(42)
                   !!! PROPS(8) beta 需要设为0
                   !!! STATEV(13*NSLPTL+28) 这个是孔隙率
                   END DO

                   TERM8=0.
                   DO K=1,3
                      TERM8=TERM8+(1.-STATEV(13*NSLPTL+28))**2.
     2                *(SLPDEFEFF(K,J)+PROPS(8)*SLPDEFNOREFF(K,J)
     3                                 *DSIGN(1.D0,FSLIP(J)))
                   !!! 这个是师兄Word的公式(42)
                   !!! PROPS(8) beta 需要设为0
                   END DO

                   WORKST(ID,J)=TERM1*TERM2*TERM6
     2                         -TERM1*TERM5*TERM7
     3                         -TERM1*TERM2*TERM4*TERM8

                   VAR1=0.
                   DO K=1,NSLPTL
                      VAR1=VAR1+SLIPGM(J,K)*STATEV(11*NSLPTL+27+K)
                   END DO
                   !!! 加上强度系数后的位错密度的和
                   WORKST(ID,J)=-TERM1*DSIGN(1.D0,FSLIP(J))
     2                /PROPS(95+(P-1)*9)
     3                *(VAR1**0.5/PROPS(101+(P-1)*9)
     4                  -2.*PROPS(98+(P-1)*9)*STATEV(11*NSLPTL+27+J))
     5                *DFDXRHO(ID,J)+WORKST(ID,J)
                   !!! PROPS(95+(P-1)*9) 柏氏矢量
                   !!! PROPS(101+(P-1)*9) 增殖系数K
                   !!! PROPS(98+(P-1)*9) 湮灭系数yc
                   IF (EQSTRESS.NE.0.) THEN
                      TERM7=0.
                      DO K=1,6
                         TERM7=TERM7+STRESS(K)*SLPDEFEFF(K,J)
                         !!! (42)式的最后一项，注意使得beta=0
                      END DO
                      
                      WORKST(ID,J)=-TERM1*TERM2*TERM4
     2                             *PROPS(14)/(PROPS(15)*(2.*PI)**0.5)
     3                             *EXP(-0.5*((STATEV(13*NSLPTL+29)
     4                                      -PROPS(16))/PROPS(15))**2.)
     5                             *TERM7/EQSTRESS
     6                             +WORKST(ID,J)
                      !!! PROPS(14) 孔洞形核系数fN
                      !!! PROPS(15) 孔洞形核系数sN
                      !!! STATEV(13*NSLPTL+29) 等效塑性应变
                      !!! PROPS(16) 孔洞形核系数eN
                      !!! EQSTRESS von-Mises应力
                   END IF

                END DO
                WORKST(ID,ID)=1.+WORKST(ID,ID)
             END DO      
          END DO
      
      ELSE
          !!! 这一部分没有变化
          !!! 本质是在求解位错随着各个滑移系
          !!! 应变的增量
          
          ID=0
          DO P=1,NSET
             DO Q=1,NSLIP(P)

                ID=ID+1
                DO J=1,NSLPTL
                   VAR1=0.
                   DO K=1,NSLPTL
                      VAR1=VAR1+SLIPGM(ID,K)*STATEV(11*NSLPTL+27+K)
                   END DO

                   DDRDDGA(ID,J)=-ABS(DGAMOD(ID))/PROPS(95+(P-1)*9)
     2                       *1./(2.*PROPS(101+(P-1)*9)
     3                              *VAR1**0.5)
     4                       *SLIPGM(ID,J)
                END DO

                DDRDDGA(ID,ID)=1.+2.*PROPS(98+(P-1)*9)/PROPS(95+(P-1)*9)
     2                           *ABS(DGAMOD(ID))+DDRDDGA(ID,ID)
             END DO
          END DO
      
  
          ID=0
          DO P=1,NSET
             DO Q=1,NSLIP(P)
             
                ID=ID+1
                DO J=1,NSLPTL
                   IF (ID.EQ.J) THEN
                       VAR1=0.
                       DO K=1,NSLPTL
                          VAR1=VAR1+SLIPGM(ID,K)*STATEV(11*NSLPTL+27+K)
                       END DO

                       DDRDDGC(ID,J)=DSIGN(1.D0,DGAMOD(ID))
     2                         /PROPS(95+(P-1)*9)
     3                         *(VAR1**0.5/PROPS(101+(P-1)*9)
     4                           -2.*PROPS(98+(P-1)*9)
     5                              *STATEV(11*NSLPTL+27+ID))
                   ELSE
                       DDRDDGC(ID,J)=0.
                   END IF
                END DO
                
             END DO
          END DO
        
        
        
          PRTSTR='DDRDDGA'
          CALL LUDCMP (DDRDDGA, NSLPTL, ND, INDX, DDCMP,PRTSTR,IDPNEWDT)
          IF (IDPNEWDT.EQ.1) THEN
             PNEWDT=0.5
             GO TO 2002
          END IF
          DO K=1,NSLPTL
             CALL LUBKSB (DDRDDGA, NSLPTL, ND, INDX, DDRDDGC(1,K))
          END DO
      

          !!! 这里的输出
          !!! 就是位错密度的增量对滑移系滑移量的增量的导数

          VAR1=0.
          DO I=1,NSLPTL
             DO J=1,3
                VAR1=VAR1+DGAMOD(I)*
     2          (SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                                  *DSIGN(1.D0,FSLIP(I)))
          !!! PROPS(8) 设置为0
             END DO
          END DO

          VAR4=0.
          !!! 这一项目与孔隙率增量对滑移系上滑移量增量的导数有关
          DO J=1,NSLPTL
             DO K=1,6
                VAR4=VAR4+STRESS(K)*DGAMOD(J)
     2          *(SLPDEFEFF(K,J)+PROPS(8)*SLPDEFNOREFF(K,J)
     3                           *DSIGN(1.D0,FSLIP(J)))
             END DO
          END DO

          DO I=1,NSLPTL
             VAR2=0.
             DO J=1,3
                VAR2=VAR2+SLPDEFEFF(J,I)
     2          +PROPS(8)*SLPDEFNOREFF(J,I)*DSIGN(1.D0,FSLIP(I))
             END DO
             
             DDFVOIDDDG(I)=(1.-STATEV(13*NSLPTL+28))**2.*VAR2
             !!! 孔隙率增量对滑移系上滑移量增量的导数
             !!! STATEV(13*NSLPTL+28) 孔隙率
             !!! VAR2 [ ]:I
             IF (EQSTRESS.NE.0.) THEN
                 VAR3=0.
                 DO J=1,6
                    VAR3=VAR3+STRESS(J)
     2              *(SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                                       *DSIGN(1.D0,FSLIP(I)))
                 END DO
                 !!! VAR3 S:[ ]
                 VAR5=EXP(-0.5*((STATEV(13*NSLPTL+29)
     2                   -PROPS(16))/PROPS(15))**2.)
                 !!! PROPS(16) 孔隙率形核系数eN
                 !!! PROPS(15) 孔隙率形核系数sN
                 !!! PROPS(14) 孔隙率形核系数fN
                 !!! STATEV(13*NSLPTL+29) 等效塑性应变
                 DDFVOIDDDG(I)=DDFVOIDDDG(I)
     2           -PROPS(14)*(STATEV(13*NSLPTL+29)-PROPS(16))
     3            /(PROPS(15)**3.*(2.*PI)**0.5)
     4            *VAR5*VAR3/EQSTRESS*VAR4/EQSTRESS
     5           +PROPS(14)/(PROPS(15)*(2.*PI)**0.5)
     6            *VAR5*VAR3/EQSTRESS

             END IF

             DDFVOIDDDG(I)=DDFVOIDDDG(I)
     2           /(1.+2.*(1.-STATEV(13*NSLPTL+28))*VAR1)
             !!! 第10页，除掉一个系数
          END DO


          TERM1=THETA*DTIME
          ID=0
          DO P=1,NSET
             DO Q=1,NSLIP(P)
                ID=ID+1

                TERM2=DFDXTAUEFF(ID)
                TERM3=DTAUEFFDTAU(ID)
                TERM4=DTAUEFFDFVOID(ID)
                TERM5=DFDXSIGEFF(ID)

                DO J=1,NSLPTL
                   TERM6=0.
                   DO K=1,6
                      TERM6=TERM6+(1.-STATEV(13*NSLPTL+28))
     2                *DDEMSDEFF0(K,ID)*(SLPDEFEFF(K,J)
     3                 +PROPS(8)*SLPDEFNOREFF(K,J)*DSIGN(1.D0,FSLIP(J)))
                      !!! 
                   END DO

                   TERM7=0.
                   DO K=1,NSLPTL
                      DO L=1,6
                         TERM7=TERM7+DDEMSDEFF0(L,ID)*(SLPDEFEFF(L,K)
     2                         +PROPS(8)*SLPDEFNOREFF(L,K)
     3                         *DSIGN(1.D0,FSLIP(K)))*DGAMOD(K)
                      END DO
                   END DO

                   TERM8=0.
                   DO K=1,6
                      TERM8=TERM8+(1.-STATEV(13*NSLPTL+28))
     2                *DDEMSDNOREFF(K,ID)*(SLPDEFEFF(K,J)
     3                 +PROPS(8)*SLPDEFNOREFF(K,J)*DSIGN(1.D0,FSLIP(J)))
                   END DO

                   TERM9=0.
                   DO K=1,NSLPTL
                      DO L=1,6
                         TERM9=TERM9+DDEMSDNOREFF(L,ID)*(SLPDEFEFF(L,K)
     2                         +PROPS(8)*SLPDEFNOREFF(L,K)
     3                         *DSIGN(1.D0,FSLIP(K)))*DGAMOD(K)
                      END DO
                   END DO

                   WORKST(ID,J)=TERM1*TERM2*TERM6
     2                         -TERM1*TERM2*DDFVOIDDDG(J)*TERM7
     3                         -TERM1*TERM5*TERM8
     4                         +TERM1*TERM5*DDFVOIDDDG(J)*TERM9
     5                         -TERM1*TERM2*TERM4*DDFVOIDDDG(J)
                   !!! 以上这些都是新增的
                   !!! F\alpha对滑移系增量的导数(t+dt)的时间下
                   VAR1=0.
                   DO K=1,NSLPTL
                      VAR1=VAR1+DFDXRHO(ID,K)*DDRDDGC(K,J)
                   END DO

                   WORKST(ID,J)=-TERM1*VAR1+WORKST(ID,J)

                END DO

                WORKST(ID,ID)=1.+WORKST(ID,ID)
             END DO
          END DO

      END IF





      PRTSTR='WORKST'
      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP,PRTSTR, IDPNEWDT)
      IF (IDPNEWDT.EQ.1) THEN
          PNEWDT=0.5
          GO TO 2002
      END IF





C-----  Increment of shear strain in a slip system: DGAMMA
      TERM1=THETA*DTIME
      DO I=1,NSLPTL
C-----EQ.38
         IF (NITRTN.EQ.0) THEN

            TERM2=DFDXTAUEFF(I)
            TERM3=DTAUEFFDTAU(I)
            !!! TEMR3好像没有用到
            TERM4=DFDXSIGEFF(I)

            VAR1=0.
            DO J=1,NDI
               VAR1=VAR1+DDEMSDEFF0(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  VAR1=VAR1+DDEMSDEFF0(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF


            VAR2=0.
            DO J=1,NDI
               VAR2=VAR2+DDEMSDNOREFF(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  VAR2=VAR2+DDEMSDNOREFF(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF
            !!! 式子42的右边一项
            DGAMMA(I)=FSLIP(I)*DTIME+TERM1*TERM2*VAR1
     2                              -TERM1*TERM4*VAR2

C-----EQ.40
         ELSE
            DGAMMA(I)=TERM1*(FSLIP(I)-FSLIP1(I))+FSLIP1(I)*DTIME
     2                -DGAMOD(I)

         END IF

      END DO

      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)


      DO I=1,NSLPTL
         DGAMMA(I)=DGAMMA(I)+DGAMOD(I)
      END DO  


C    calculate DRHO in each slip system    
       IF (NITRTN.EQ.0) THEN
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)
                ID=ID+1
                VAR1=0.
                DO K=1,NSLPTL
                   VAR1=VAR1+SLIPGM(ID,K)
     2                       *(STATEV(11*NSLPTL+27+K)-DRHOOD(K))
                END DO

                DRHOSP(ID)=1./PROPS(95+(I-1)*9)
     2                     *(VAR1**0.5/PROPS(101+(I-1)*9)
     3                -2.*PROPS(98+(I-1)*9)
     4                   *(STATEV(11*NSLPTL+27+ID)-DRHOOD(ID)))
     5                  *ABS(DGAMMA(ID))
                !!! 位错密度增量
                !!! PROPS(95+(I-1)*9) 柏氏矢量
                !!! PROPS(101+(I-1)*9) 位错增殖系数K
                !!! PROPS(98+(I-1)*9) 位错湮灭系数K
             END DO
          END DO

       ELSE

C---- 4th order Runge-Kutta method  
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)
                ID=ID+1
                VAR1=0.
                DO K=1,NSLPTL
                   VAR1=VAR1+SLIPGM(ID,K)
     2                       *(STATEV(11*NSLPTL+27+K)-DRHOOD(K))
                END DO

                DRHORK45(ID,1)=1./PROPS(95+(I-1)*9)
     2            *(VAR1**0.5/PROPS(101+(I-1)*9)
     3              -2.*PROPS(98+(I-1)*9)*(STATEV(11*NSLPTL+27+ID)
     4                                              -DRHOOD(ID)))
     5            *ABS(DGAMMA(ID))
             END DO
          END DO

          
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)
                ID=ID+1
                VAR1=0.
                DO K=1,NSLPTL
                   VAR1=VAR1+SLIPGM(ID,K)
     2                *(STATEV(11*NSLPTL+27+K)-DRHOOD(K)
     3                                +DRHORK45(K,1)/3.)
                END DO

                DRHORK45(ID,2)=1./PROPS(95+(I-1)*9)
     2          *(VAR1**0.5/PROPS(101+(I-1)*9)
     3            -2.*PROPS(98+(I-1)*9)
     4             *(STATEV(11*NSLPTL+27+ID)-DRHOOD(ID)
     5                              +DRHORK45(ID,1)/3.))
     6            *ABS(DGAMMA(ID))
             END DO
          END DO

          
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)
                ID=ID+1
                VAR1=0.
                DO K=1,NSLPTL
                   VAR1=VAR1+SLIPGM(ID,K)
     2            *(STATEV(11*NSLPTL+27+K)-DRHOOD(K)-DRHORK45(K,1)/3.
     3                  +DRHORK45(K,2))
                END DO
                DRHORK45(ID,3)=1./PROPS(95+(I-1)*9)
     2            *(VAR1**0.5/PROPS(101+(I-1)*9)
     3              -2.*PROPS(98+(I-1)*9)
     4            *(STATEV(11*NSLPTL+27+ID)-DRHOOD(ID)-DRHORK45(ID,1)/3.
     5                                                 +DRHORK45(ID,2)))
     6            *ABS(DGAMMA(ID))
             END DO
          END DO

  
          ID=0
          DO I=1,NSET
             DO J=1,NSLIP(I)
                ID=ID+1
                VAR1=0.
                DO K=1,NSLPTL
                   VAR1=VAR1+SLIPGM(ID,K)
     2                  *(STATEV(11*NSLPTL+27+K)-DRHOOD(K)
     3                                     +DRHORK45(K,3))
                END DO
                DRHORK45(ID,4)=1./PROPS(95+(I-1)*9)
     2             *(VAR1**0.5/PROPS(101+(I-1)*9)
     3               -2.*PROPS(98+(I-1)*9)
     4                *(STATEV(11*NSLPTL+27+ID)-DRHOOD(ID)
     5                                    +DRHORK45(ID,3)))
     6             *ABS(DGAMMA(ID))               
             END DO
          END DO

          DO I=1,NSLPTL
             DRHOSP(I)=(DRHORK45(I,1)+3.*DRHORK45(I,2)
     2                  +3.*DRHORK45(I,3)+DRHORK45(I,4))/8.
          END DO
          
      END IF


C-----  Update the shear strain in a slip system: STATEV(NSLPTL+1) - 
C       STATEV(2*NSLPTL)
C-----  Update the cumulative shear strain in a slip system: 
C       STATEV(9*NSLPTL+1) - STATEV(10*NSLPTL)
      
      DO I=1,NSLPTL
         STATEV(NSLPTL+I)=STATEV(NSLPTL+I)+DGAMMA(I)-DGAMOD(I)
         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))
     2                      -ABS(DGAMOD(I))
      END DO


C-----  Update the current dislocation in a slip system: STATEV(11*NSLPTL+27+1) - 
C     STATEV(11*NSLPTL)
      
      STATEV(11*NSLPTL+20)=0.
      DO I=1,NSLPTL
         STATEV(11*NSLPTL+27+I)=STATEV(11*NSLPTL+27+I)
     2                          +DRHOSP(I)-DRHOOD(I)
         !!! 各个滑移系上的位错密度
         STATEV(11*NSLPTL+20)=STATEV(11*NSLPTL+20)
     3                        +STATEV(11*NSLPTL+27+I)
         !!! 总的位错密度
         IF (STATEV(11*NSLPTL+27+I).LE.0.) THEN
             STATEV(11*NSLPTL+27+I)=0.
         END IF
         !!! 位错密度不能小于零
      END DO


      
      

      
      
C-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-(1.-STATEV(13*NSLPTL+28))
     2      *(SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                       *DSIGN(1.D0,FSLIP(I)))*DGAMMA(I)
         END DO
            !!! 减去塑性应变的增量
            !!! 得到弹性变形张量
            !!! 控制PROPS(8)是个零
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-(1.-STATEV(13*NSLPTL+28))
     2      *(SLPDEFEFF(J+3,I)+PROPS(8)*SLPDEFNOREFF(J+3,I)
     3                         *DSIGN(1.D0,FSLIP(I)))*DGAMMA(I)
         END DO
      END DO



C---- Update elastic, plastic and total strain components

      IF (NITRTN.EQ.0) THEN
          CALL ROTSIG(STATEV(11*NSLPTL+1),DROT,ELASTRAN,2,NDI,NSHR)
      END IF

      DO I=1,NTENS
         ELASTRAN(I)=ELASTRAN(I)+DELATS(I)-DELATSOD(I)
      END DO


      IF (NITRTN.EQ.0) THEN
          CALL ROTSIG(STATEV(11*NSLPTL+7),DROT,PLASTRAN,2,NDI,NSHR)
      END IF

      DO I=1,NTENS
         DPLATS(I)=DSTRAN(I)-DELATS(I)
         PLASTRAN(I)=PLASTRAN(I)+DPLATS(I)-DPLATSOD(I)
      END DO

      DO I=1,NTENS
         STATEV(11*NSLPTL+12+I)=ELASTRAN(I)+PLASTRAN(I)
      END DO

C---- Update equivalent plastic strain
      STATEV(11*NSLPTL+19)=((PLASTRAN(1)**2.+PLASTRAN(2)**2.
     2                       +PLASTRAN(3)**2.+0.5*PLASTRAN(4)**2.
     3                       +0.5*PLASTRAN(5)**2.+0.5*PLASTRAN(6)**2.)
     4                       *2/3)**0.5
      !!! STATEV(11*NSLPTL+19) 等效塑性应变
C---- Update cumulative equivalent plastic strain
      DEQPLATS=DPLATS(1)**2.+DPLATS(2)**2.+DPLATS(3)**2.
     2         +0.5*DPLATS(4)**2.+0.5*DPLATS(5)**2.+0.5*DPLATS(6)**2.
     
      DEQPLATS=(DEQPLATS*2./3.)**0.5
      
      STATEV(NSTATV-5)=STATEV(NSTATV-5)+DEQPLATS-DEQPLATSOD
      !!! STATEV(NSTATV-5) 等效塑性应变增量





C---- Update tr(Dp)
      TRDP=DPLATS(1)+DPLATS(2)+DPLATS(3)

C---- Update microscopic cumulative equivalent plastic strain
C---- and volume fraction of void

      DFVOID=(1.-STATEV(13*NSLPTL+28))*TRDP


C      IF (EQSTRESS.NE.0.) THEN
C
C          VAR1=0.
C          DO I=1,NSLPTL
C             DO J=1,6
C                VAR1=VAR1+STRESS(J)*DGAMMA(I)
C     2               *(SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
C     3                                     *DSIGN(1.D0,FSLIP(I)))
C             END DO
C          END DO
C      
C          DMEQPLATS=VAR1/EQSTRESS



          DO J=1,6
             DPLATSM(J)=0.
          END DO
          
          DO J=1,3
             DO I=1,NSLPTL
                DPLATSM(J)=DPLATSM(J)+
     2          (SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                           *DSIGN(1.D0,FSLIP(I)))*DGAMMA(I)
             END DO
          END DO
          
          DO J=1,3
             DO I=1,NSLPTL
                DPLATSM(J+3)=DPLATSM(J+3)+
     2          (SLPDEFEFF(J+3,I)+PROPS(8)*SLPDEFNOREFF(J+3,I)
     3                             *DSIGN(1.D0,FSLIP(I)))*DGAMMA(I)
             END DO
          END DO
          
          DMEQPLATS=DPLATSM(1)**2.+DPLATSM(2)**2.+DPLATSM(3)**2.
     2         +0.5*DPLATSM(4)**2.+0.5*DPLATSM(5)**2.+0.5*DPLATSM(6)**2.


          DMEQPLATS=(2./3.*DMEQPLATS)**0.5

          IF (ISNAN(DMEQPLATS)) THEN
              DMEQPLATS=10.**10.
          END IF

          DFVOID1=PROPS(14)/(PROPS(15)*(2.*PI)**0.5)
     2             *EXP(-0.5*((STATEV(13*NSLPTL+29)
     3                  -DMEQPLATSOD-PROPS(16))
     4                  /PROPS(15))**2.)*DMEQPLATS
          
          STATEV(13*NSLPTL+29)=STATEV(13*NSLPTL+29)
     2                         +DMEQPLATS-DMEQPLATSOD

          DFVOID2=PROPS(14)/(PROPS(15)*(2.*PI)**0.5)
     2             *EXP(-0.5*((STATEV(13*NSLPTL+29)-PROPS(16))
     3                  /PROPS(15))**2.)*DMEQPLATS

          DFVOID=DFVOID+(DFVOID1+DFVOID2)/2.

          STATEV(13*NSLPTL+28)=STATEV(13*NSLPTL+28)+DFVOID-DFVOIDOD
          !!! STATEV(13*NSLPTL+28) 孔隙率
          !!! DFVOID 孔隙率增量当前值
          !!! DFVOIDOD 孔隙率增量过去值
          IF (STATEV(13*(NSLPTL)+28).GE.1.) THEN
              STATEV(13*(NSLPTL)+28)=1.
          ELSE IF((STATEV(13*(NSLPTL)+28).LE.0.)) THEN
              STATEV(13*(NSLPTL)+28)=0.
          END IF


     
C      ELSE
C
CC          DMEQPLATS=0.
C
C          DMEQPLATS=DEQPLATS/(1-STATEV(13*NSLPTL+28))
C
C          DFVOID=DFVOID
C
C          STATEV(13*NSLPTL+28)=STATEV(13*NSLPTL+28)+DFVOID-DFVOIDOD
C
C          IF (STATEV(13*(NSLPTL)+28).GE.1.) THEN
C              STATEV(13*(NSLPTL)+28)=1.
C          ELSE IF((STATEV(13*(NSLPTL)+28).LE.0.)) THEN
C              STATEV(13*(NSLPTL)+28)=0.
C          END IF
C
C          STATEV(13*NSLPTL+29)=STATEV(13*NSLPTL+29)
C     2                         +DMEQPLATS-DMEQPLATSOD
C
C      END IF
      







C-----  Increment of deformation gradient associated with lattice 
C     stretching in the current state, i.e. the velocity gradient 
C     (associated with lattice stretching) times the increment of time:
C     DVGRAD (only needed for finite rotation)
C     De+W-Wp=De-We=1/2*(Le+Transpose(Le))-1/2*(Le-Transpose(Le))=Transpose(Le)
C     Eqs. (8), (9) and (10)
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.
                  ELSE
                     TERM1=-1.
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPNEFF(IJ2,K)*
     3                                (1.-STATEV(13*NSLPTL+28))
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPNEFF(IJ2,K)*
     3                                (1.-STATEV(13*NSLPTL+28))
                  END DO
               END IF

            END DO
         END DO

      END IF



C-----  Increment of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.
         END DO
      ELSE
         DO I=1,NTENS
C----    original (HUANG 1991, Hill and Rice 1972)
           DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF



      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-(1.-STATEV(13*NSLPTL+28))
     2        *(DDEMSDEFF(I,J)+PROPS(8)*DDEMSDNOREFF(I,J)
     3                        *DSIGN(1.D0,FSLIP(J)))*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR

            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)
     2           -(DDEMSDEFF(I+3,J)+PROPS(8)*DDEMSDNOREFF(I+3,J)
     3                           *DSIGN(1.D0,FSLIP(J)))*DGAMMA(J)
            END DO

         END DO
      END IF



C-----  Update the stress: STRESS
      DO I=1,NTENS
         STRESS(I)=STRESS(I)+DSTRES(I)-DSOLD(I)
      END DO



C-----  Increment of normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.
               DSPDIR(I,J)=0.

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO


C-----  Update the normal to a slip plane and a slip direction (only 
C     needed for finite rotation)
C
         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)-DSPNRO(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)-DSPDRO(I,J)
            END DO
         END DO

      END IF

C-----Slip deformation tensor: SLPDEF (Schmid factors)
      DO J=1,NSLPTL
         SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
         SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
         SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
         SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
         SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
         SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         
         SLPDEFNOR(1,J)=SLPNOR(1,J)*SLPNOR(1,J)
         SLPDEFNOR(2,J)=SLPNOR(2,J)*SLPNOR(2,J)
         SLPDEFNOR(3,J)=SLPNOR(3,J)*SLPNOR(3,J)
         SLPDEFNOR(4,J)=SLPNOR(1,J)*SLPNOR(2,J)
     2                 +SLPNOR(2,J)*SLPNOR(1,J)
         SLPDEFNOR(5,J)=SLPNOR(1,J)*SLPNOR(3,J)
     2                 +SLPNOR(3,J)*SLPNOR(1,J)
         SLPDEFNOR(6,J)=SLPNOR(2,J)*SLPNOR(3,J)
     2                 +SLPNOR(3,J)*SLPNOR(2,J)
      END DO



C-----Update value of the resolved shear stress in slip systems
      DO I=1,NSLPTL
         TERM1=0.

         DO J=1,NTENS
            IF (J.LE.NDI) THEN
               TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
            ELSE
               TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
            END IF
         END DO

         STATEV(2*NSLPTL+I)=TERM1
      END DO



C     calculate Mises stress, hydrostatic pressure, deviatoric stress
      HYDROSTRESS=(STRESS(1)+STRESS(2)+STRESS(3))/3.
      DEVSIG(1,1)=STRESS(1)-HYDROSTRESS
      DEVSIG(2,2)=STRESS(2)-HYDROSTRESS
      DEVSIG(3,3)=STRESS(3)-HYDROSTRESS
      DEVSIG(1,2)=STRESS(4)
      DEVSIG(1,3)=STRESS(5)
      DEVSIG(2,3)=STRESS(6)
      DEVSIG(2,1)=STRESS(4)
      DEVSIG(3,1)=STRESS(5)
      DEVSIG(3,2)=STRESS(6)
      EQSTRESS=0.
      DO I=1,3
         DO J=1,3
            EQSTRESS=EQSTRESS+DEVSIG(I,J)*DEVSIG(I,J)
         END DO
      END DO

      EQSTRESS=(3./2.*EQSTRESS)**0.5

C---- Update effective void fraction
      IF (STATEV(13*(NSLPTL)+28).GE.(1./PROPS(12))) THEN
          FVOIDEFF=1./PROPS(12)
      ELSE IF (STATEV(13*(NSLPTL)+28).LE.0.) THEN
          FVOIDEFF=0.
      ELSE
          FVOIDEFF=STATEV(13*(NSLPTL)+28)
      END IF

C---- Update effective resolved shear and normal stress
      DO I=1,NSLPTL
         IDTAUEFF(I)=0
      END DO

      IF (EQSTRESS.EQ.0.) THEN
         DO I=1,NSLPTL
            TAUSPEFF(I)=0.
         END DO
      ELSE
         DO I=1,NSLPTL
            IF (STATEV(2*NSLPTL+I).EQ.0.) THEN
                TAUSPEFF(I)=0.
            ELSE
                ID=-1
9003            CONTINUE
                ID=ID+1
                IF (ID.EQ.0) THEN
                   TAUSPEFF(I)=abs(STATEV(2*NSLPTL+I))
                ELSE
                   VAR1=(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     2                 +2.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
     3                    *DSIGN(1.0,STATEV(2*NSLPTL+I))
     4                    /(TAUSPEFF(I))**2.
     2                 +2./45.*PROPS(11)*FVOIDEFF
     3                        *(EQSTRESS/TAUSPEFF(I))**2.
     4                 +2.*PROPS(12)*FVOIDEFF
     5                   *(1.+0.5*PROPS(13)**2.*(
     6                   3./10.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**2.+
     7                   3./20.*(HYDROSTRESS/TAUSPEFF(I))**2.))
     7                 -1.-(PROPS(12)*FVOIDEFF)**2.

                   VAR3=4.*1./TAUSPEFF(I)
                   VAR3=VAR3*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
                   VAR3=VAR3*DSIGN(1.0,STATEV(2*NSLPTL+I))
                   VAR3=VAR3/TAUSPEFF(I)**2.

                   VAR2=-2.*1./TAUSPEFF(I)
     2                     *(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     3                  -VAR3
     3                  -4./45.*PROPS(11)*FVOIDEFF
     4                    *1./TAUSPEFF(I)*(EQSTRESS/TAUSPEFF(I))**2.
     5                  -2.*PROPS(12)*FVOIDEFF*0.5*PROPS(13)**2.*
     6                   (3./5.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**3.
     7                   +3./10.*HYDROSTRESS**2./TAUSPEFF(I)**3.)

                   TAUSPEFF(I)=TAUSPEFF(I)-VAR1/VAR2
                END IF

                RESIDU=(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     2                 +2.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)
     3                    *DSIGN(1.0,STATEV(2*NSLPTL+I))
     4                    /(TAUSPEFF(I))**2.
     2                 +2./45.*PROPS(11)*FVOIDEFF
     3                        *(EQSTRESS/TAUSPEFF(I))**2.
     4                 +2.*PROPS(12)*FVOIDEFF
     5                   *(1.+0.5*PROPS(13)**2.*(
     6                   3./10.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**2.+
     7                   3./20.*(HYDROSTRESS/TAUSPEFF(I))**2.))
     7                 -1.-(PROPS(12)*FVOIDEFF)**2.

                IF (ID.EQ.0) THEN
                    GO TO 9003
                END IF

                IF (ISNAN(RESIDU).OR.(ABS(RESIDU).GT.(10.**10.))) THEN
                   TAUSPEFF(I)=STATEV(2*NSLPTL+I)
                   IDTAUEFF(I)=1
                   GO TO 9004
                ELSE IF (ABS(RESIDU).GT.(10.**-6.)) THEN
                   GO TO 9003
                ELSE IF (ABS(RESIDU).LE.(10.**-6.)) THEN
                   IDTAUEFF(I)=2
                   GO TO 9004
                END IF
            END IF

            TAUSPEFF(I) = TAUSPEFF(I)*DSIGN(1.0,STATEV(2*NSLPTL+I))

9004        CONTINUE
         END DO

      END IF

      SLPDEFEFF33=0.

      IF (EQSTRESS.EQ.0.) THEN
          DO I=1,NSLPTL
             DO J=1,3
                DO K=1,3
                   SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                END DO
             END DO
          END DO
      ELSE
          DO I=1,NSLPTL
             IF (IDTAUEFF(I).EQ.0) THEN
                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                    END DO
                 END DO
             ELSE IF (IDTAUEFF(I).EQ.1) THEN
                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=SLPDIR(J,I)*SLPNOR(K,I)
                    END DO
                 END DO
             ELSE
                 VAR1=-2.*1./TAUSPEFF(I)
     2                     *(STATEV(2*NSLPTL+I)/TAUSPEFF(I))**2.
     3                  -4.*FVOIDEFF*PROPS(17)*STATEV(2*NSLPTL+I)*
     4                     *DSIGN(1.0,STATEV(2*NSLPTL+I))
     5                     /TAUSPEFF(I)**3.
     3                  -4./45.*PROPS(11)*FVOIDEFF
     4                    *1./TAUSPEFF(I)*(EQSTRESS/TAUSPEFF(I))**2.
     5                  -2.*PROPS(12)*FVOIDEFF*0.5*PROPS(12)**2.*
     6                   (3./5.*PROPS(17)*HYDROSTRESS/TAUSPEFF(I)**3.
     7                   +3./10.*HYDROSTRESS**2./TAUSPEFF(I)**3.)
                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=2.*STATEV(2*NSLPTL+I)
     2                                      /TAUSPEFF(I)**2.
     3                                     *SLPDIR(J,I)*SLPNOR(K,I)
     4                                   +2.*FVOIDEFF*PROPS(17)
     5                                      *DSIGN(1.0,STATEV(2*NSLPTL+I))
     6                                      /TAUSPEFF(I)**2.
     7                                     *SLPDIR(J,I)*SLPNOR(K,I)
     4                                   +2./15.*PROPS(11)*FVOIDEFF
     5                                  /TAUSPEFF(I)**2.*DEVSIG(J,K)
                    END DO
                    SLPDEFEFF33(J,J,I)=SLPDEFEFF33(J,J,I)
     5                   +2.*PROPS(12)*FVOIDEFF*0.5*PROPS(13)**2.
     6                    /10.*(PROPS(17)+HYDROSTRESS)
     7                    /TAUSPEFF(I)**2.
                 END DO

                 DO J=1,3
                    DO K=1,3
                       SLPDEFEFF33(J,K,I)=-SLPDEFEFF33(J,K,I)/VAR1
                    END DO
                 END DO

             END IF
          END DO
      END IF

      SLPDEFNOREFF33=0.
      DO I=1,NSLPTL
         DO J=1,3
            VAR1=SLPDEFEFF33(J,1,I)**2+SLPDEFEFF33(J,2,I)**2
     2           +SLPDEFEFF33(J,3,I)**2
            IF (VAR1.NE.0.) THEN
               DO K=1,3
                  DO L=1,3
                     SLPDEFNOREFF33(K,L,I)=
     2               SLPDEFEFF33(J,K,I)*SLPDEFEFF33(J,L,I)/VAR1
                  END DO
               END DO
               GO TO 7002
            END IF
         END DO
7002     CONTINUE
      END DO


      SLPDEFEFF=0.
      DO I=1,NSLPTL
         SLPDEFEFF(1,I)=SLPDEFEFF33(1,1,I)
         SLPDEFEFF(2,I)=SLPDEFEFF33(2,2,I)
         SLPDEFEFF(3,I)=SLPDEFEFF33(3,3,I)
         SLPDEFEFF(4,I)=SLPDEFEFF33(1,2,I)+SLPDEFEFF33(2,1,I)
         SLPDEFEFF(5,I)=SLPDEFEFF33(1,3,I)+SLPDEFEFF33(3,1,I)
         SLPDEFEFF(6,I)=SLPDEFEFF33(2,3,I)+SLPDEFEFF33(3,2,I)
      END DO


      SLPDEFNOREFF=0.
      DO I=1,NSLPTL
         SLPDEFNOREFF(1,I)=SLPDEFNOREFF33(1,1,I)
         SLPDEFNOREFF(2,I)=SLPDEFNOREFF33(2,2,I)
         SLPDEFNOREFF(3,I)=SLPDEFNOREFF33(3,3,I)
         SLPDEFNOREFF(4,I)=SLPDEFNOREFF33(1,2,I)
     2                    +SLPDEFNOREFF33(2,1,I)
         SLPDEFNOREFF(5,I)=SLPDEFNOREFF33(1,3,I)
     2                    +SLPDEFNOREFF33(3,1,I)
         SLPDEFNOREFF(6,I)=SLPDEFNOREFF33(2,3,I)
     2                    +SLPDEFNOREFF33(3,2,I)
      END DO



C---- Update effective resolved normal stress
      DO I=1,NSLPTL
         TERM1=0.
      
         DO J=1,NTENS
            IF (J.LE.NDI) THEN
               TERM1=TERM1+SLPDEFNOREFF(J,I)*STRESS(J)
            ELSE
               TERM1=TERM1+SLPDEFNOREFF(J-NDI+3,I)*STRESS(J)
            END IF
         END DO
      
         STATEV(10*NSLPTL+I)=-TERM1+PROPS(10)
      END DO





C-----  Update the current strength in a slip system: STATEV(1) - 
C     STATEV(NSLPTL)
      ID=0
      DO I=1,NSET
         DO J=1,NSLIP(I)
            ID=ID+1
            VAR1=0.
            DO K=1,NSLPTL
               VAR1=VAR1+SLIPDM(ID,K)*STATEV(11*NSLPTL+27+K)
            END DO
            STATEV(ID)=(PROPS(93+(I-1)*9)
     2                  +PROPS(92)*SHEARMD
     3                   *PROPS(95+(I-1)*9)*(VAR1**0.5))
     4                 *(1.+PROPS(9)*STATEV(10*NSLPTL+ID))
         END DO
      END DO







C-----  Derivative of shear strain increment in a slip system w.r.t. 
C     strain increment: DDGDDE
C
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
            TERM2=DFDXTAUEFF(J)
            TERM3=DTAUEFFDTAU(J)
            TERM4=DFDXSIGEFF(J)
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM1*TERM2*DDEMSDEFF0(I,J)
            ELSE
               DDGDDE(J,I)=TERM1*TERM2*DDEMSDEFF0(I-NDI+3,J)
            END IF
            
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=DDGDDE(J,I)-TERM1*TERM4*DDEMSDNOREFF(I,J)
            ELSE
               DDGDDE(J,I)=DDGDDE(J,I)
     2                     -TERM1*TERM4*DDEMSDNOREFF(I-NDI+3,J)
            END IF
            
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO



       VAR1=0.
       DO I=1,NSLPTL
          DO J=1,3
             VAR1=VAR1+DGAMMA(I)*
     2       (SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                               *DSIGN(1.D0,FSLIP(I)))
          END DO
       END DO

       VAR4=0.
       DO J=1,NSLPTL
          DO K=1,6
             VAR4=VAR4+STRESS(K)*DGAMMA(J)
     2       *(SLPDEFEFF(K,J)+PROPS(8)*SLPDEFNOREFF(K,J)
     3                             *DSIGN(1.D0,FSLIP(J)))
          END DO
       END DO


       DO M=1,NTENS
          DDFVOIDDDE(M)=0.
          DO I=1,NSLPTL
             VAR2=0.
             DO J=1,3
                VAR2=VAR2+SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     2                                           *DSIGN(1.D0,FSLIP(I))
             END DO

             VAR2=(1.-STATEV(13*NSLPTL+28))**2.*VAR2
             
             IF (EQSTRESS.NE.0.) THEN
                 VAR3=0.
                 DO J=1,6
                    VAR3=VAR3+STRESS(J)
     2              *(SLPDEFEFF(J,I)+PROPS(8)*SLPDEFNOREFF(J,I)
     3                                     *DSIGN(1.D0,FSLIP(J)))
                 END DO

                 VAR5=EXP(-0.5*((STATEV(13*NSLPTL+29)
     2                   -PROPS(16))/PROPS(15))**2.)

                 VAR2=VAR2
     2           -PROPS(14)*(STATEV(13*NSLPTL+29)-PROPS(16))
     3            /(PROPS(15)**3.*(2.*PI)**0.5)
     4            *VAR5*VAR3/EQSTRESS*VAR4/EQSTRESS
     5           +PROPS(14)/(PROPS(15)*(2.*PI)**0.5)
     6            *VAR5*VAR3/EQSTRESS
             END IF

             VAR2=VAR2*DDGDDE(I,M)
             DDFVOIDDDE(M)=DDFVOIDDDE(M)
     2              +VAR2/(1.+2.*(1.-STATEV(13*NSLPTL+28))*VAR1)
          END DO
       END DO









C-----  Derivative of stress increment w.r.t. strain increment, i.e. 
C     Jacobian matrix
C
C-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF



C-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-(1.-STATEV(13*NSLPTL+28))
     2         *(DDEMSDEFF(I,K)+PROPS(8)*DDEMSDNOREFF(I,K)
     3                          *DSIGN(1.D0,FSLIP(K)))*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2           (1.-STATEV(13*NSLPTL+28))*(DDEMSDEFF(I+3,K)
     3               +PROPS(8)*DDEMSDNOREFF(I+3,K)*DSIGN(1.D0,FSLIP(K)))
     4               *DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2            (1.-STATEV(13*NSLPTL+28))
     3            *(DDEMSDEFF(I,K)+PROPS(8)*DDEMSDNOREFF(I,K)
     4              *DSIGN(1.D0,FSLIP(K)))*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2            (1.-STATEV(13*NSLPTL+28))
     3            *(DDEMSDEFF(J+3,K)+PROPS(8)*DDEMSDNOREFF(J+3,K)
     4              *DSIGN(1.D0,FSLIP(K)))*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF



      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL

               DDSDDE(I,J)=DDSDDE(I,J)
     2         +(DDEMSDEFF(I,K)+PROPS(8)
     3           *DDEMSDNOREFF(I,K)*DSIGN(1.D0,FSLIP(K)))
     4           *DGAMMA(K)*DDFVOIDDDE(J)


            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL

                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)+
     2            (DDEMSDEFF(I+3,K)+PROPS(8)
     3             *DDEMSDNOREFF(I+3,K)*DSIGN(1.D0,FSLIP(K)))
     4             *DGAMMA(K)*DDFVOIDDDE(J+NDI)
                   !!! DDFVOIDDDE
                   !!! 孔隙率的增量对应变张量增量的导数
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)
     2            +(DDEMSDEFF(I,K)+PROPS(8)
     3              *DDEMSDNOREFF(I,K)*DSIGN(1.D0,FSLIP(K)))
     4              *DDFVOIDDDE(J+NDI)*DGAMMA(K)

                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)
     2            +(DDEMSDEFF(J+3,K)+PROPS(8)
     3              *DDEMSDNOREFF(J+3,K)*DSIGN(1.D0,FSLIP(K)))
     4              *DDFVOIDDDE(I)*DGAMMA(K)
               END DO
            END DO

         END DO
      END IF





      IF (ITRATN.NE.0) THEN
         DO J=1,NTENS
            DO I=1,NTENS
               DDSDDE(I,J)=DDSDDE(I,J)/(1.+DEV)
            END DO
         END DO
      END IF

      
    


C-----  Iteration ?
      IF (ITRATN.NE.0) THEN

C-----  Save solutions (without iteration):
C            Shear strain-rate in a slip system FSLIP1
C            Current strength in a slip system GSLP1
C            Shear strain in a slip system GAMMA1
C            Resolved shear stress in a slip system TAUSP1
C            Normal to a slip plane SPNOR1
C            Slip direction SPDIR1
C            Stress STRES1
C            Jacobian matrix DDSDE1
C
         IF (NITRTN.EQ.0) THEN

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               FSLIP1(J)=FSLIP(J)
               GSLP1(J)=STATEV(J)
               GAMMA1(J)=STATEV(NSLPTL+J)
               TAUSP1(J)=STATEV(2*NSLPTL+J)
               SGAMMA1(J)=STATEV(9*NSLPTL+J)
               SIGSPEFF1(J)=STATEV(10*NSLPTL+J)
               RHOSP1(J)=STATEV(11*NSLPTL+27+J)
               TAUSPEFF1(J)=TAUSPEFF(J)
               
               DO I=1,3
                  IDNOR=IDNOR+1
                  SPNOR1(I,J)=STATEV(IDNOR)

                  IDDIR=IDDIR+1
                  SPDIR1(I,J)=STATEV(IDDIR)
               END DO
            END DO

            DO J=1,NTENS
               STRES1(J)=STRESS(J)
               ELASTRAN1(J)=ELASTRAN(J)
               PLASTRAN1(J)=PLASTRAN(J)
               TOTALSTRAN1(J)=STATEV(11*NSLPTL+12+J)
               DO I=1,NTENS
                  DDSDE1(I,J)=DDSDDE(I,J)
               END DO
            END DO
            
            TRDP1=TRDP
            
            EQPLATS1=STATEV(11*NSLPTL+19)
            CEQPLATS1=STATEV(NSTATV-5)
            FVOID1=STATEV(13*NSLPTL+28)
            EQPLATSM1=STATEV(13*NSLPTL+29)
            
         END IF

C-----  Increments of stress DSOLD, and solution dependent state 
C     variables DGAMOD, DTAUOD, DRHOOD, DSPNRO, DSPDRO (for the next 
C     iteration)
C
         DO I=1,NTENS
            DSOLD(I)=DSTRES(I)
         END DO

         DO J=1,NSLPTL
            DGAMOD(J)=DGAMMA(J)
            DTAUOD(J)=DTAUSP(J)
            DRHOOD(J)=DRHOSP(J)
            DO I=1,3
               DSPNRO(I,J)=DSPNOR(I,J)
               DSPDRO(I,J)=DSPDIR(I,J)
            END DO
         END DO

         DO I=1,NTENS
            DELATSOD(I)=DELATS(I)
            DPLATSOD(I)=DPLATS(I)
         END DO

         DEQPLATSOD=DEQPLATS
         DFVOIDOD=DFVOID
         DMEQPLATSOD=DMEQPLATS

C-----  Check if the iteration solution converges
         IDBACK=0
         ID=0
         DO I=1,NSET
            DO J=1,NSLIP(I)
               ID=ID+1
               VAR1=0.
               DO K=1,NSLPTL
                  VAR1=VAR1+SLIPDM(ID,K)*STATEV(11*NSLPTL+27+K)
               END DO
               VAR1=(PROPS(93+(I-1)*9)
     2               +PROPS(92)*SHEARMD*PROPS(95+(I-1)*9)
     3               *VAR1**0.5)*(1.+PROPS(9)*STATEV(10*NSLPTL+ID))
               VAR1=PROPS(99+(I-1)*9)
     2              *(ABS(TAUSPEFF(ID)/VAR1))
     3                **PROPS(100+(I-1)*9)
     4              *DSIGN(1.D0,TAUSPEFF(ID))
               RESIDU=THETA*DTIME*VAR1
     2                +DTIME*(1.0-THETA)*FSLIP1(ID)-DGAMMA(ID)
               
               
               IF (ABS(RESIDU).GT.GAMERR) IDBACK=1
            END DO
         END DO

         IF ((IDBACK.NE.0).AND.(NITRTN.LT.ITRMAX)) THEN

            GO TO 2001



         ELSE IF (NITRTN.GE.ITRMAX) THEN
C-----  Solution not converge within maximum number of iteration (the 
C     solution without iteration will be used)
C
            DO J=1,NTENS
               STRESS(J)=STRES1(J)
               DO I=1,NTENS
                  DDSDDE(I,J)=DDSDE1(I,J)
               END DO
            END DO

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               STATEV(J)=GSLP1(J)
               STATEV(NSLPTL+J)=GAMMA1(J)
               STATEV(2*NSLPTL+J)=TAUSP1(J)
               STATEV(9*NSLPTL+J)=SGAMMA1(J)
               STATEV(10*NSLPTL+J)=SIGSPEFF1(J)
               STATEV(11*NSLPTL+27+J)=RHOSP1(J)
               TAUSPEFF(J)=TAUSPEFF1(J)

               DO I=1,3
                  IDNOR=IDNOR+1
                  STATEV(IDNOR)=SPNOR1(I,J)

                  IDDIR=IDDIR+1
                  STATEV(IDDIR)=SPDIR1(I,J)
               END DO
            END DO

            DO I=1,NTENS
               ELASTRAN(I)=ELASTRAN1(I)
               PLASTRAN(I)=PLASTRAN1(I)
               STATEV(11*NSLPTL+12+I)=TOTALSTRAN1(I)
            END DO

            TRDP=TRDP1

            STATEV(11*NSLPTL+19)=EQPLATS1
            STATEV(NSTATV-5)=CEQPLATS1
            STATEV(13*NSLPTL+28)=FVOID1
            STATEV(13*NSLPTL+29)=EQPLATSM1

         END IF

      END IF

C     Calculate the determinant of deformation gradient 
      STATEV(NSTATV-4)=DFGRD1(1,1)*DFGRD1(2,2)*DFGRD1(3,3)
     1                 +DFGRD1(1,2)*DFGRD1(2,3)*DFGRD1(3,1)
     2                 +DFGRD1(1,3)*DFGRD1(3,2)*DFGRD1(2,1)
     3                 -DFGRD1(1,2)*DFGRD1(2,1)*DFGRD1(3,3) 
     4                 -DFGRD1(1,3)*DFGRD1(3,1)*DFGRD1(2,2)
     5                 -DFGRD1(2,3)*DFGRD1(3,2)*DFGRD1(1,1)

C---- Elastic and plastic strain components
      DO I=1,NTENS
         STATEV(11*NSLPTL+I)=ELASTRAN(I)
         STATEV(11*NSLPTL+6+I)=PLASTRAN(I)
      END DO
      
C---- Upadate effective resolved shear stress
      DO I=1,NSLPTL
         STATEV(12*NSLPTL+27+I)=TAUSPEFF(I)
      END DO

C---- Upadate equivalent strain
      STATEV(11*(NSLPTL)+27)=
     1    (2./3*( STATEV(11*(NSLPTL)+12+1)**2.
     2           +STATEV(11*(NSLPTL)+12+2)**2.
     3           +STATEV(11*(NSLPTL)+12+3)**2.
     4           +0.5*STATEV(11*(NSLPTL)+12+4)**2.
     5           +0.5*STATEV(11*(NSLPTL)+12+5)**2.
     6           +0.5*STATEV(11*(NSLPTL)+12+6)**2.))**0.5

      STATEV(NSTATV-6)=NITRTN


C---- Plastic volume increase
      STATEV(NSTATV-8)=STATEV(NSTATV-8)+STATEV(NSTATV-4)*TRDP

CC---- update effective void fraction
      IF (STATEV(13*(NSLPTL)+28).GE.(1./PROPS(12))) THEN
          STATEV(NSTATV-7)=1./PROPS(12)
      ELSE IF (STATEV(13*(NSLPTL)+28).LE.0.) THEN
          STATEV(NSTATV-7)=0.
      ELSE
          STATEV(NSTATV-7)=STATEV(13*(NSLPTL)+28)
      END IF

2002  CONTINUE
      RETURN
      END
C----------------------------------------------------------------------





C----------------------------------------------------------------------
      SUBROUTINE ROTATION (PROP, ROTATE)

C-----  This subroutine calculates the rotation matrix, i.e. the 
C     direction cosines of cubic crystal [100], [010] and [001] 
C     directions in global system

C-----  The rotation matrix is stored in the array ROTATE.

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3) 
      CHARACTER*80 PRTSTR
C-----  Subroutines:
C
C       CROSS  -- cross product of two vectors
C
C       LUDCMP -- LU decomposition
C
C       LUBKSB -- linear equation solver based on LU decomposition 
C                 method (must call LUDCMP first)


C-----  PROP -- constants characterizing the crystal orientation 
C               (INPUT)
C
C            PROP(1) - PROP(3) -- direction of the first vector in 
C                                 local cubic crystal system
C            PROP(4) - PROP(6) -- direction of the first vector in 
C                                 global system
C
C            PROP(9) - PROP(11)-- direction of the second vector in 
C                                 local cubic crystal system
C            PROP(12)- PROP(14)-- direction of the second vector in 
C                                 global system
C
C-----  ROTATE -- rotation matrix (OUTPUT):
C
C            ROTATE(i,1) -- direction cosines of direction [1 0 0] in 
C                           local cubic crystal system
C            ROTATE(i,2) -- direction cosines of direction [0 1 0] in 
C                           local cubic crystal system
C            ROTATE(i,3) -- direction cosines of direction [0 0 1] in 
C                           local cubic crystal system

C-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)

C-----  LU decomposition of TERM1
      PRTSTR='ROTATION'
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP,PRTSTR, IDPNEWDT)

C-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.
            ELSE
               TERM2(I,J)=0.
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

C-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)

C-----  Check: the angle between first and second vector in local and 
C     global systems must be the same.  The relative difference must be
C     less than 0.1%.
C
      IF (ABS(ANGLE1/ANGLE2-1.).GT.0.001) THEN 
         WRITE (6,*) 
     2      '***ERROR - angles between two vectors are not the same'
         STOP
      END IF

C-----  rotation matrix: ROTATE
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END
C----------------------------------------------------------------------



C----------------------------------------------------------------------
           SUBROUTINE CROSS (A, B, C, ANGLE)

C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              STOP
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*) 
     2           '***ERROR - first and second vectors are parallel'
               STOP
            END IF

           RETURN
           END


C----------------------------------------------------------------------



C----------------------------------------------------------------------
           SUBROUTINE CROSS1 (A, B, C, ANGLE)

C-----  (1) normalize vectors A and B to unit vectors
C       (2) store A, B and A*B (cross product) in C

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           DO I=1,3
              C(I,1)=A(I)
              C(I,2)=B(I)
           END DO

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)


           RETURN
           END
C----------------------------------------------------------------------





C----------------------------------------------------------------------
      SUBROUTINE SLIPSYS (ISPDIR, ISPNOR, NSLIP, SLPDIR, SLPNOR, 
     2                    ROTATE)

C-----  This subroutine generates all slip systems in the same set for 
C     a CUBIC crystal.  For other crystals (e.g., HCP, Tetragonal, 
C     Orthotropic, ...), it has to be modified to include the effect of
C     crystal aspect ratio.

C-----  Denote s as a slip direction and m as normal to a slip plane.  
C     In a cubic crystal, (s,-m), (-s,m) and (-s,-m) are NOT considered
C     independent of (s,m).

C-----  Subroutines:  LINE1 and LINE

C-----  Variables:
C
C     ISPDIR -- a typical slip direction in this set of slip systems 
C               (integer)  (INPUT)
C     ISPNOR -- a typical normal to slip plane in this set of slip 
C               systems (integer)  (INPUT)
C     NSLIP  -- number of independent slip systems in this set 
C               (OUTPUT)
C     SLPDIR -- unit vectors of all slip directions  (OUTPUT)
C     SLPNOR -- unit normals to all slip planes  (OUTPUT)
C     ROTATE -- rotation matrix (INPUT)
C          ROTATE(i,1) -- direction cosines of [100] in global system
C          ROTATE(i,2) -- direction cosines of [010] in global system
C          ROTATE(i,3) -- direction cosines of [001] in global system
C
C     NSPDIR -- number of all possible slip directions in this set
C     NSPNOR -- number of all possible slip planes in this set
C     IWKDIR -- all possible slip directions (integer)
C     IWKNOR -- all possible slip planes (integer)


C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ISPDIR(3), ISPNOR(3), SLPDIR(3,50), SLPNOR(3,50), 
     *          ROTATE(3,3), IWKDIR(3,24), IWKNOR(3,24), TERM(3)

      NSLIP=0
      NSPDIR=0
      NSPNOR=0

C-----  Generating all possible slip directions in this set
C
C       Denote the slip direction by [lmn].  I1 is the minimum of the 
C     absolute value of l, m and n, I3 is the maximum and I2 is the 
C     mode, e.g. (1 -3 2), I1=1, I2=2 and I3=3.  I1<=I2<=I3.

      I1=MIN(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I3=MAX(IABS(ISPDIR(1)),IABS(ISPDIR(2)),IABS(ISPDIR(3)))
      I2=IABS(ISPDIR(1))+IABS(ISPDIR(2))+IABS(ISPDIR(3))-I1-I3

              RMODIR=SQRT(FLOAT(I1*I1+I2*I2+I3*I3))

C     I1=I2=I3=0
      IF (I3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip direction is [000]'
         STOP

C     I1=I2=0, I3>0   ---   [001] type
      ELSE IF (I2.EQ.0) THEN
         NSPDIR=3
         DO J=1,3
            DO I=1,3
               IWKDIR(I,J)=0
               IF (I.EQ.J) IWKDIR(I,J)=I3
            END DO
         END DO

C     I1=0, I3>=I2>0
      ELSE IF (I1.EQ.0) THEN

C        I1=0, I3=I2>0   ---   [011] type
         IF (I2.EQ.I3) THEN
            NSPDIR=6
            DO J=1,6
               DO I=1,3
                  IWKDIR(I,J)=I2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKDIR(I,J)=0
                  IWKDIR(1,6)=-I2
                  IWKDIR(2,4)=-I2
                  IWKDIR(3,5)=-I2
               END DO
            END DO

C        I1=0, I3>I2>0   ---   [012] type
         ELSE
            NSPDIR=12
            CALL LINE1 (I2, I3, IWKDIR(1,1), 1)
            CALL LINE1 (I3, I2, IWKDIR(1,3), 1)
            CALL LINE1 (I2, I3, IWKDIR(1,5), 2)
            CALL LINE1 (I3, I2, IWKDIR(1,7), 2)
            CALL LINE1 (I2, I3, IWKDIR(1,9), 3)
            CALL LINE1 (I3, I2, IWKDIR(1,11), 3)

         END IF

C     I1=I2=I3>0   ---   [111] type
      ELSE IF (I1.EQ.I3) THEN
         NSPDIR=4
         CALL LINE (I1, I1, I1, IWKDIR)

C     I3>I2=I1>0   ---   [112] type
      ELSE IF (I1.EQ.I2) THEN
         NSPDIR=12
         CALL LINE (I1, I1, I3, IWKDIR(1,1))
         CALL LINE (I1, I3, I1, IWKDIR(1,5))
         CALL LINE (I3, I1, I1, IWKDIR(1,9))

C     I3=I2>I1>0   ---   [122] type
      ELSE IF (I2.EQ.I3) THEN
         NSPDIR=12
         CALL LINE (I1, I2, I2, IWKDIR(1,1))
         CALL LINE (I2, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I2, I1, IWKDIR(1,9))

C     I3>I2>I1>0   ---   [123] type
      ELSE
         NSPDIR=24
         CALL LINE (I1, I2, I3, IWKDIR(1,1))
         CALL LINE (I3, I1, I2, IWKDIR(1,5))
         CALL LINE (I2, I3, I1, IWKDIR(1,9))
         CALL LINE (I1, I3, I2, IWKDIR(1,13))
         CALL LINE (I2, I1, I3, IWKDIR(1,17))
         CALL LINE (I3, I2, I1, IWKDIR(1,21))

      END IF

C-----  Generating all possible slip planes in this set
C
C       Denote the normal to slip plane by (pqr).  J1 is the minimum of
C     the absolute value of p, q and r, J3 is the maximum and J2 is the
C     mode, e.g. (1 -2 1), J1=1, J2=1 and J3=2.  J1<=J2<=J3.

      J1=MIN(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J3=MAX(IABS(ISPNOR(1)),IABS(ISPNOR(2)),IABS(ISPNOR(3)))
      J2=IABS(ISPNOR(1))+IABS(ISPNOR(2))+IABS(ISPNOR(3))-J1-J3

      RMONOR=SQRT(FLOAT(J1*J1+J2*J2+J3*J3))

      IF (J3.EQ.0) THEN 
         WRITE (6,*) '***ERROR - slip plane is [000]'
         STOP

C     (001) type
      ELSE IF (J2.EQ.0) THEN
         NSPNOR=3
         DO J=1,3
            DO I=1,3
               IWKNOR(I,J)=0
               IF (I.EQ.J) IWKNOR(I,J)=J3
            END DO
         END DO

      ELSE IF (J1.EQ.0) THEN

C     (011) type
         IF (J2.EQ.J3) THEN
            NSPNOR=6
            DO J=1,6
               DO I=1,3
                  IWKNOR(I,J)=J2
                  IF (I.EQ.J.OR.J-I.EQ.3) IWKNOR(I,J)=0
                  IWKNOR(1,6)=-J2
                  IWKNOR(2,4)=-J2
                  IWKNOR(3,5)=-J2
               END DO
            END DO

C     (012) type
         ELSE
            NSPNOR=12
            CALL LINE1 (J2, J3, IWKNOR(1,1), 1)
            CALL LINE1 (J3, J2, IWKNOR(1,3), 1)
            CALL LINE1 (J2, J3, IWKNOR(1,5), 2)
            CALL LINE1 (J3, J2, IWKNOR(1,7), 2)
            CALL LINE1 (J2, J3, IWKNOR(1,9), 3)
            CALL LINE1 (J3, J2, IWKNOR(1,11), 3)

         END IF

C     (111) type
      ELSE IF (J1.EQ.J3) THEN
         NSPNOR=4
         CALL LINE (J1, J1, J1, IWKNOR)

C     (112) type
      ELSE IF (J1.EQ.J2) THEN
         NSPNOR=12
         CALL LINE (J1, J1, J3, IWKNOR(1,1))
         CALL LINE (J1, J3, J1, IWKNOR(1,5))
         CALL LINE (J3, J1, J1, IWKNOR(1,9))

C     (122) type
      ELSE IF (J2.EQ.J3) THEN
         NSPNOR=12
         CALL LINE (J1, J2, J2, IWKNOR(1,1))
         CALL LINE (J2, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J2, J1, IWKNOR(1,9))

C     (123) type
      ELSE
         NSPNOR=24
         CALL LINE (J1, J2, J3, IWKNOR(1,1))
         CALL LINE (J3, J1, J2, IWKNOR(1,5))
         CALL LINE (J2, J3, J1, IWKNOR(1,9))
         CALL LINE (J1, J3, J2, IWKNOR(1,13))
         CALL LINE (J2, J1, J3, IWKNOR(1,17))
         CALL LINE (J3, J2, J1, IWKNOR(1,21))

      END IF

C-----  Generating all slip systems in this set
C
C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in local cubic crystal system
C
C      WRITE (6,*) '          '
C      WRITE (6,*) ' #          Slip plane          Slip direction'

      DO J=1,NSPNOR
         DO I=1,NSPDIR

            IDOT=0
            DO K=1,3
               IDOT=IDOT+IWKDIR(K,I)*IWKNOR(K,J)
            END DO

            IF (IDOT.EQ.0) THEN
               NSLIP=NSLIP+1
               DO K=1,3
                  SLPDIR(K,NSLIP)=IWKDIR(K,I)/RMODIR
                  SLPNOR(K,NSLIP)=IWKNOR(K,J)/RMONOR
               END DO

C               WRITE (6,10) NSLIP, 
C     2                      (IWKNOR(K,J),K=1,3), (IWKDIR(K,I),K=1,3)

            END IF

         END DO
      END DO
10    FORMAT(1X,I2,9X,'(',3(1X,I2),1X,')',10X,'[',3(1X,I2),1X,']')

C      WRITE (6,*) 'Number of slip systems in this set = ',NSLIP
C      WRITE (6,*) '          '

      IF (NSLIP.EQ.0) THEN
         WRITE (6,*) 
     *      'There is no slip direction normal to the slip planes!'
         STOP

      ELSE

C-----  Unit vectors in slip directions: SLPDIR, and unit normals to 
C     slip planes: SLPNOR in global system
C
         DO J=1,NSLIP
            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPDIR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPDIR(I,J)=TERM(I)
            END DO

            DO I=1,3
               TERM(I)=0.
               DO K=1,3
                  TERM(I)=TERM(I)+ROTATE(I,K)*SLPNOR(K,J)
               END DO
            END DO
            DO I=1,3
               SLPNOR(I,J)=TERM(I)
            END DO
         END DO

      END IF

      RETURN
      END
C----------------------------------------------------------------------



C----------------------------------------------------------------------
           SUBROUTINE LINE (I1, I2, I3, IARRAY)

C-----  Generating all possible slip directions <lmn> (or slip planes 
C     {lmn}) for a cubic crystal, where l,m,n are not zeros.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,4)

           DO J=1,4
              IARRAY(1,J)=I1
              IARRAY(2,J)=I2
              IARRAY(3,J)=I3
           END DO

           DO I=1,3
              DO J=1,4
                 IF (J.EQ.I+1) IARRAY(I,J)=-IARRAY(I,J)
              END DO
           END DO

           RETURN
           END
C----------------------------------------------------------------------



C----------------------------------------------------------------------
           SUBROUTINE LINE1 (J1, J2, IARRAY, ID)

C-----  Generating all possible slip directions <0mn> (or slip planes 
C     {0mn}) for a cubic crystal, where m,n are not zeros and m does 
C     not equal n.

C-----  Use single precision on cray
C
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION IARRAY(3,2)

           IARRAY(ID,1)=0
           IARRAY(ID,2)=0

           ID1=ID+1
           IF (ID1.GT.3) ID1=ID1-3
           IARRAY(ID1,1)=J1
           IARRAY(ID1,2)=J1

           ID2=ID+2
           IF (ID2.GT.3) ID2=ID2-3
           IARRAY(ID2,1)=J2
           IARRAY(ID2,2)=-J2
  
           RETURN
           END
C----------------------------------------------------------------------




C----------------------------------------------------------------------
      SUBROUTINE LUDCMP (A, N, NP, INDX, D, PRTSTR, IDPNEWDT)

C-----  LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200, TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)
      DIMENSION FSLIP(NP), DFDXTAU(NP), DFDXRHO(NP),DRHOSP(NP)
      CHARACTER*80 PRTSTR

      D=1.
      DO I=1,N
         AAMAX=0.

         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         END DO

         IF (AAMAX.EQ.0.) THEN
C             write(6,*)  'Singular matrix.'
             write(6,*)  PRTSTR
             write(6,*)  A
C             STOP 
             IDPNEWDT=1
             GO TO 8001
         ELSE
             IDPNEWDT=0
         END IF
 
         VV(I)=1./AAMAX
      END DO

      DO J=1,N
         DO I=1,J-1
            SUM=A(I,J)

            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
         END DO
         AAMAX=0.

         DO I=J,N
            SUM=A(I,J)

            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            END DO

            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            END IF
         END DO

         IF (J.NE.IMAX) THEN
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            END DO

            D=-D
            VV(IMAX)=VV(J)
         END IF

         INDX(J)=IMAX
         IF (A(J,J).EQ.0.) A(J,J)=TINY
         IF (J.NE.N) THEN
            DUM=1./A(J,J)
            DO I=J+1,N
               A(I,J)=A(I,J)*DUM
            END DO
         END IF

      END DO

8001  CONTINUE

      RETURN
      END
C----------------------------------------------------------------------
      
      
      
C----------------------------------------------------------------------
      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

C-----  Linear equation solver based on LU decomposition

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      II=0
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF
         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END

C----------------------------------------------------------------------



C----------------------------------------------------------------------      
      SUBROUTINE SLIPSYSCOE (CMT,ND,NSLPTL,PROP)

C-----  Initialization for the interaction coefficient matrix

C-----  Use single precision on cray
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMT(ND,ND), PROP(6)
 
      DO I=1,NSLPTL
         CMT(I,I)=PROP(1)
      END DO
        
      CMT(1,2)=PROP(2)
      CMT(1,3)=PROP(2)
      CMT(1,4)=PROP(5)
      CMT(1,5)=PROP(5)
      CMT(1,6)=PROP(4)
      CMT(1,7)=PROP(3)
      CMT(1,8)=PROP(6)
      CMT(1,9)=PROP(5)
      CMT(1,10)=PROP(3)
      CMT(1,11)=PROP(6)
      CMT(1,12)=PROP(5)
      CMT(2,3)=PROP(2)
      CMT(2,4)=PROP(3)
      CMT(2,5)=PROP(6)
      CMT(2,6)=PROP(5)
      CMT(2,7)=PROP(5)
      CMT(2,8)=PROP(5)
      CMT(2,9)=PROP(4)
      CMT(2,10)=PROP(6)
      CMT(2,11)=PROP(3)
      CMT(2,12)=PROP(5)
      CMT(3,4)=PROP(6)
      CMT(3,5)=PROP(3)
      CMT(3,6)=PROP(5)
      CMT(3,7)=PROP(6)
      CMT(3,8)=PROP(3)
      CMT(3,9)=PROP(5)
      CMT(3,10)=PROP(5)
      CMT(3,11)=PROP(5)
      CMT(3,12)=PROP(4)
      CMT(4,5)=PROP(2)
      CMT(4,6)=PROP(2)
      CMT(4,7)=PROP(6) 
      CMT(4,8)=PROP(5)
      CMT(4,9)=PROP(3)
      CMT(4,10)=PROP(5)
      CMT(4,11)=PROP(4)
      CMT(4,12)=PROP(5)
      CMT(5,6)=PROP(2)
      CMT(5,7)=PROP(5)
      CMT(5,8)=PROP(4)
      CMT(5,9)=PROP(5)
      CMT(5,10)=PROP(6)
      CMT(5,11)=PROP(5)
      CMT(5,12)=PROP(3)
      CMT(6,7)=PROP(3)
      CMT(6,8)=PROP(5)
      CMT(6,9)=PROP(6)
      CMT(6,10)=PROP(3)
      CMT(6,11)=PROP(5)
      CMT(6,12)=PROP(6)
      CMT(7,8)=PROP(6)
      CMT(7,9)=PROP(2)
      CMT(7,10)=PROP(4)
      CMT(7,11)=PROP(5)
      CMT(7,12)=PROP(5)
      CMT(8,9)=PROP(2)
      CMT(8,10)=PROP(5)
      CMT(8,11)=PROP(6)
      CMT(8,12)=PROP(3)
      CMT(9,10)=PROP(5)
      CMT(9,11)=PROP(3)
      CMT(9,12)=PROP(6)
      CMT(10,11)=PROP(2)
      CMT(10,12)=PROP(2)
      CMT(11,12)=PROP(2)
    
      DO I=2,NSLPTL
         DO J=1,(I-1)
             CMT(I,J)=CMT(J,I)
         END DO
      END DO
      
      RETURN
      END
C----------------------------------------------------------------------