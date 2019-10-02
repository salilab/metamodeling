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
glu_measurement = [92, 102, 130, 155, 169, 173, 173, 171, 167, 160, 150, 140, 132, 124, 117, 110, 107, 102, 100, 100, 99, 99, 98, 98, 97, 97, 96, 96, 96, 96, 95, 95, 95, 92, 92, 90, 90, 90, 90, 90, 90, 90, 90]

print(len(glu_measurement))

glucose = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement]

#print(glucose)

glucose1=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 11)[0:10])

print(glucose1)
print(len(glucose1))

#----------write out the measurement----------
f1=open("meal_exp1_normal2.dat","w")

for x in zip(time, glucose1):
    f1.write("{0}\t{1}\n".format(*x))
