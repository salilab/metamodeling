#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
from math import *
import csv
import operator
import random

def intakefunct(x):
    glucose = 10/(1+np.exp(6-0.2*x))
    dg = 2*np.exp(x/5+6)/(np.exp(x/5)+np.exp(6))**2
    return glucose, dg

# 420 min in total, one measurement per 20 min, 21 timepoints
time = np.arange(420*1)
#print(np.linspace(0, 60, 61))
glucose=[]
dg=[]
glucose.append(0)
for i in range(11,70):
    glucose.append(intakefunct(i)[0])
    dg.append(intakefunct(i+1)[0]-intakefunct(i)[0])
for i in range(70,431):
    glucose.append(intakefunct(59)[0])
    dg.append(0)
print(glucose)
print(dg)

#----------write out the measurement----------
f1=open("meal_exp1_normal2_Gin_dt1_sigmoid.dat","w")

for x in zip(time, dg, glucose):
    f1.write("{0}\t{1}\t{2}\n".format(*x))
    
