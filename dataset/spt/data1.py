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
glucose = [25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5]

glucose1=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 10))

#print(glucose1)
print(len(glucose1))

glucose2 = random.sample(range(70, 200), 21)
#print(glucose2)

# plasma insulin concentration after a meal - pmol/l
insulin = [70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170]

insulin1=[]
for i in range(0,len(insulin)-1):
    insulin1.extend(np.linspace(insulin[i], insulin[i+1], 10))
    
insulin2 = random.sample(range(20, 400), 21)

#----------write out the measurement----------
f1=open("spt_obs1_avr.dat","w")
f2=open("spt_obs2_random.dat","w")

for x in zip(time, glucose1, insulin1):
    f1.write("{0}\t{1}\t{2}\n".format(*x))

for x in zip(time, glucose2, insulin2):
    f2.write("{0}\t{1}\t{2}\n".format(*x))
