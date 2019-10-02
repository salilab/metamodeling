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
# 240 min in total, one measurement per 12 min, 20 timepoints
time = np.arange(140*3)

# plasma concentration of GLP1, pmol/L
glp1_normal_measurement = [i*1E-9 for i in [4.8, 6.5, 11.8, 15, 16, 17, 16.8, 17.2, 17.3, 16.8, 16, 14.3, 13.5, 12.8, 11.8, 11, 10.2, 9.5, 9, 8.6, 8.4]]
glp1_t2d_measurement = [i*1E-9 for i in [6.5, 7.2, 10.5, 13.4, 14.2, 14, 13, 12, 11.5, 11.1, 
11, 10.2, 9.8, 9.5, 9.3, 9.2, 8.8, 8.5, 8.4, 8.4, 8.4]]

glp1_normal=[]
glp1_t2d=[]
for i in range(0,len(glp1_normal_measurement)-1):
    glp1_normal.extend(np.linspace(glp1_normal_measurement[i], glp1_normal_measurement[i+1], 5)[0:4])

for i in range(80, 140):
    glp1_normal.append(glp1_normal_measurement[-1])
    
print(glp1_normal)
for i in range(0,len(glp1_t2d_measurement)-1):    
    glp1_t2d.extend(np.linspace(glp1_t2d_measurement[i], glp1_t2d_measurement[i+1], 5)[0:4])
    
for i in range(80, 140):
    glp1_t2d.append(glp1_t2d_measurement[-1])
    
print(glp1_normal)
print(len(glp1_normal))
print(glp1_t2d)
print(len(glp1_t2d))
#----------write out the measurement----------
f1=open("glp1_exp1_normal_t2d.dat","w")

for x in zip(time, glp1_normal, glp1_t2d):
    f1.write("{0}\t{1}\n".format(*x))
    
