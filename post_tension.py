#!/usr/bin/python
#-- coding:utf8 --
from abaqus import getInput,getInputs
from odbAccess import *
from abaqusConstants import *
from textRepr import *
import os
import csv

Jobname = 'Job-2'

odbFile = Jobname+'.odb'
csvFile = Jobname+'.csv'

odb = openOdb(path = odbFile)
frame_length = len(odb.steps.values()[-1].frames)

fixSet = odb.rootAssembly.instances['TESS-1']

RF33 = [0]*frame_length
LE33 = [0]*frame_length

for j in range(0,frame_length):
    field_RF = odb.steps.values()[-1].frames[j].fieldOutputs['S']
    subField_RF = field_RF.getSubset(region = fixSet)
    sum_RF = 0
    num = 0
    for val_1 in subField_RF.values:
        sum_RF = sum_RF + val_1
        num = num+1
    RF = sum_RF
    RF33[j] =  sum_RF.data[2]/num/1.e6

for k in range(0,frame_length):
    field_U = odb.steps.values()[-1].frames[k].fieldOutputs['LE']
    subField_U = field_U.getSubset(region = fixSet)
    U3 = 0
    num = 0
    for val_2 in subField_U.values:
        U3 = U3+val_2
        num = num+1
    LE33[k] = U3.data[2]/num;

temp = []

f = open(csvFile, "w")
f.write("Strain,")
f.write("Stress/MPa,")
f.write('\n')
for i in range(0,frame_length):
    f.write(str(LE33[i]))
    f.write(',')
    f.write(str(RF33[i]))
    f.write('\n')
f.close()

### OUPUT ### 