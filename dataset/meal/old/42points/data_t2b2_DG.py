#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(42*10)

# plasma glucose concentration after a meal - mg/dl
glu_measurement = [165, 170, 189, 218, 245, 270, 290, 302, 318, 328, 332, 338, 339, 338, 334, 328, 320, 310, 300, 289, 276, 266, 250, 241, 230, 220, 213, 209, 202, 199, 195, 192, 188, 181, 178, 171, 166, 162, 157, 152, 149, 147, 142]

print(len(glu_measurement))

glucose = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement]

#print(glucose)

glucose1=[]
DG=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 11)[0:10])

for i in range(0,len(glucose1)-1):
    tmp = glucose1[i+1]-glucose1[i]
    if tmp <=0:
        tmp = 0
    DG.append(tmp)

DG.append(0)
print(glucose1)
print(len(glucose1))

#----------write out the measurement----------
f1=open("meal_exp1_t2b2.dat","w")
f2=open("meal_exp1_t2b2_DG.dat","w")

for x in zip(time, glucose1):
    f1.write("{0}\t{1}\n".format(*x))
    
for x in zip(time, DG):
    f2.write("{0}\t{1}\n".format(*x))
