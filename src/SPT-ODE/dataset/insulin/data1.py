#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(51*10)

# plasma insulin concentration after a meal - pmol/l
insulin=[]
for i in range(0,51):
    j = i*10.0 
    insulin.append(j)

print(len(insulin))

#----------write out the measurement----------
f1=open("insulin_exp1.dat","w")

for x in zip(time, insulin):
    f1.write("{0}\t{1}\n".format(*x))
