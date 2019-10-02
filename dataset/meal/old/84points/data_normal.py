#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(21*10)

# plasma glucose concentration after a meal - mg/dl
glu_measurement = [92, 140, 178, 175, 160, 142, 127, 112, 102, 98, 95, 95, 95, 95, 95, 95, 94, 94, 93, 92, 91]


glucose = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement]

print(glucose)

glucose1=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 12)[0:11])

#print(glucose1)
print(len(glucose1))
glucose2 = random.sample(range(70, 200), 21)
#print(glucose2)

# plasma insulin concentration after a meal - pmol/l
insulin = [25, 230, 300, 326, 360, 315, 220, 138, 80, 58, 42, 38, 35, 34, 30, 30, 28, 26, 24, 22, 21]

insulin1=[]
for i in range(0,len(insulin)-1):
    insulin1.extend(np.linspace(insulin[i], insulin[i+1], 12)[0:11])

#print(glucose1)
print(len(insulin1))

insulin2 = random.sample(range(20, 400), 21)

#----------write out the measurement----------
f1=open("meal_exp1_normal.dat","w")
f2=open("meal_exp2_random.dat","w")

for x in zip(time, glucose1, insulin1):
    f1.write("{0}\t{1}\t{2}\n".format(*x))

for x in zip(time, glucose2, insulin2):
    f2.write("{0}\t{1}\t{2}\n".format(*x))
