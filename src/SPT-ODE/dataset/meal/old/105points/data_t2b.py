#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# 420 min in total, one measurement per 2 min, 210 timepoints
time = np.arange(21*10)

# plasma glucose concentration after a meal - mg/dl
glu_measurement = [165, 185, 250, 292, 320, 336, 338, 330, 312, 285, 268, 240, 220, 204, 198, 190, 180, 170, 158, 150, 140]

glucose = [x*0.001 /(180 *0.001 *0.1) for x in glu_measurement]

print(glucose)

glucose1=[]
for i in range(0,len(glucose)-1):
    glucose1.extend(np.linspace(glucose[i], glucose[i+1], 12)[0:11])

#print(glucose1)
print(len(glucose1))

# plasma insulin concentration after a meal - pmol/l
insulin = [52, 140, 240, 222, 240, 260, 290, 340, 370, 385, 370, 338, 290, 240, 200, 170, 146, 120, 104, 85, 70]

insulin1=[]
for i in range(0,len(insulin)-1):
    insulin1.extend(np.linspace(insulin[i], insulin[i+1], 12)[0:11])

#print(glucose1)
print(len(insulin1))

#insulin2 = random.sample(range(20, 400), 21)

#----------write out the measurement----------
f1=open("meal_exp1_t2b.dat","w")
#f2=open("meal_exp2_random2.dat","w")

for x in zip(time, glucose1, insulin1):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
