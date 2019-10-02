#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# define the function

def network_camp(camp0, time, dt):
    Y=[None]*len(time)
    Y[0] = 1
    for i in range(0, len(time)-1):
        Y[i+1] = Y[i] + (0.00005 + 3 - 0.5*Y[i]**2)* dt
    return Y

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
time = timeslice= list(range(0, 50,1))
print(time)
# plasma glucose concentration after a meal - mg/dl

camp=network_camp(0.2, time, 1)

#----------write out the measurement----------
f1=open("network_camp.dat","w")

for x in zip(time, camp):
    f1.write("{0}\t{1}\n".format(*x))
    
