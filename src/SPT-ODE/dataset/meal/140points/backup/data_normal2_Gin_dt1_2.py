#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

intake = gaussian(np.linspace(0, 90, 90), 45, 5)*6
print(intake)
# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(420*1)

# plasma glucose concentration after a meal - mg/dl, eat glucose for 60 min, 21 time points
glucose = []
dg=[]
glucose.extend(np.arange(0,6,0.1))
print(glucose)

for i in range(0,len(glucose)):
    dg.append(0.1)
for i in range(len(glucose),420):    
    glucose.append(0)
    dg.append(0)
    
print(dg)
print(glucose)

#----------write out the measurement----------
f1=open("meal_exp1_normal2_Gin_dt1_2.dat","w")

for x in zip(time, dg, glucose):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
    
