#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

### Running mean/Moving average
def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result
# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(21*10)

# plasma glucose concentration after a meal - mg/dl
glu_measurement = [92, 102, 130, 155, 169, 173, 173, 171, 167, 160, 150, 140, 132, 124, 117, 110, 107, 102, 100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 96, 96, 95, 95, 95, 92, 92, 90, 90, 90, 90, 90, 90, 90, 90]

print(len(glu_measurement))

glucose = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement]

#print(glucose)

glucose1=[]
DG=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 6)[0:5])

for i in range(0,len(glucose1)-1):
    tmp = glucose1[i+1]-glucose1[i]
    if tmp <=0:
        tmp = 0
    DG.append(tmp)

DG.append(0)
DG1=running_mean(DG, 3)
DG1.append(0)
print(glucose1)
print(len(glucose1))

#----------write out the measurement----------
f1=open("meal_exp1_normal2.dat","w")
f2=open("meal_exp1_normal2_DG.dat","w")

for x in zip(time, glucose1):
    f1.write("{0}\t{1}\n".format(*x))
    
for x in zip(time, DG1):
    f2.write("{0}\t{1}\n".format(*x))
