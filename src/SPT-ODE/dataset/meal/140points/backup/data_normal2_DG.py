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
time = np.arange(140*3)

# plasma glucose concentration after a meal - mg/dl
glu_measurement1 = [92, 95, 105, 126, 140, 152, 164, 169, 172, 174, 174, 172, 171, 170, 166, 160, 156, 150, 147, 140, 134]
glu_measurement2 = [134, 124, 117,  109, 103, 101, 100, 99, 97, 97, 95, 95, 95, 95, 95, 94, 94, 92, 92, 90, 90, 90, 90, 90, 89, 89]

insulin_measurement1 = [25, 230, 300, 326, 360, 315, 220, 138, 80, 58, 42, 38, 35, 34, 30, 30, 28, 26, 24, 22, 21]
insulin_measurement2 = [52, 140, 240, 222, 240, 260, 290, 340, 370, 385, 370, 338, 290, 240, 200, 170, 146, 120, 104, 85, 70]

glucose1 = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement1]
glucose2 = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement2]

print(glucose1)
print(glucose2)

glucose=[]
DG=[]
insulin=[]
for i in range(0,len(glucose1)-1):
    glucose.extend(np.linspace(glucose1[i], glucose1[i+1], 3)[0:2])
for i in range(0,len(glucose2)-1):    
    glucose.extend(np.linspace(glucose2[i], glucose2[i+1], 5)[0:4])
    
for i in range(0,len(insulin_measurement1)-1):
    insulin.extend(np.linspace(insulin_measurement1[i], insulin_measurement1[i+1], 8)[0:7])

for i in range(0,len(glucose)-1):
    tmp = glucose[i+1]-glucose[i]
    if tmp <=0:
        tmp = 0
    DG.append(tmp)

DG.append(0)
DG1=running_mean(DG, 3)
print(glucose)
print(len(glucose))

#----------write out the measurement----------
f1=open("meal_exp1_normal2.dat","w")
f2=open("meal_exp1_normal2_DG.dat","w")

for x in zip(time, glucose, insulin):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
    
for x in zip(time, DG1):
    f2.write("{0}\t{1}\n".format(*x))
