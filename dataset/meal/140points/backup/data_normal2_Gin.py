#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(140*3)
print(np.arange(21*0.5))
# plasma glucose concentration after a meal - mg/dl, eat glucose for 60 min, 21 time points
glucose = []
dg=[]
glucose.extend(np.arange(21*0.5))

for i in range(0,len(glucose)):
    dg.append(0.2)
for i in range(len(glucose),140):    
    glucose.append(0)
    dg.append(0)
    
print(dg)
print(glucose)

#----------write out the measurement----------
f1=open("meal_exp1_normal2_Gin.dat","w")

for x in zip(time, dg, glucose):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
    
