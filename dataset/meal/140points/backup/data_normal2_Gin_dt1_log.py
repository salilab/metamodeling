#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

def intakefunct(time):
    glucose = np.log(1+5*time**2)*10/np.log(1+5*60**2)
    dg = 2*time/(1+5*time**2)*10/np.log(1+5*60**2)
    return glucose, dg

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(420*1)
#print(np.linspace(0, 60, 61))
glucose=[]
dg=[]
glucose.append(0)
for i in range(1,60):
    glucose.append(intakefunct(i)[0])
    dg.append(intakefunct(i+1)[0]-intakefunct(i)[0])
for i in range(60,421):
    glucose.append(0)
    dg.append(0)
print(glucose)
print(dg)

#----------write out the measurement----------
f1=open("meal_exp1_normal2_Gin_dt1_log.dat","w")

for x in zip(time, dg, glucose):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
    
