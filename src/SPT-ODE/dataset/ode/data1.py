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
glucose = [85, 85, 90, 110, 135, 160, 175, 180, 175, 160, 150, 125, 110, 95, 90, 90, 90, 90, 90, 90, 90]

glucose1=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 10))

#print(glucose1)
print(len(glucose1))
glucose2 = random.sample(range(70, 200), 21)
#print(glucose2)

# plasma insulin concentration after a meal - pmol/l
insulin = [25, 45, 55, 110, 180, 350, 305, 310, 360, 300, 330, 250, 170, 100, 80, 75, 70, 65, 65, 50, 45]

insulin1=[]
for i in range(0,len(insulin)-1):
    insulin1.extend(np.linspace(insulin[i], insulin[i+1], 10))

#print(glucose1)
print(len(insulin1))

insulin2 = random.sample(range(20, 400), 21)

#----------write out the measurement----------
f1=open("ode_exp1_avr.dat","w")
f2=open("ode_exp2_random.dat","w")

for x in zip(time, glucose1, insulin1):
    f1.write("{0}\t{1}\t{2}\n".format(*x))

for x in zip(time, glucose2, insulin2):
    f2.write("{0}\t{1}\t{2}\n".format(*x))
