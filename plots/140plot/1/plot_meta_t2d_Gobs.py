#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc
import matplotlib
import matplotlib.ticker as mtick
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rcParams
from pylab import *
import pylab as pylab
import matplotlib.patches as patches
from matplotlib.pyplot import *
import csv
from scipy.signal import savgol_filter

# Read data file

with open("../../bnet/meta_S_t2d.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()
    
meal_G=[]
meal_Gstd=[]
meal_S=[]
meal_Sstd=[]
spt_S=[]
spt_Sstd=[]
network_S=[]
network_Sstd=[]
meta_S=[]
meta_Sstd=[]
meta_S2=[]
meta_Sstd2=[]

for i in range(0, len(data)):
    meal_G.append(float(data[i][0]))
    meal_Gstd.append(float(data[i][2]))
    meal_S.append(float(data[i][3]))
    meal_Sstd.append(float(data[i][5]))
    spt_S.append(float(data[i][6]))
    spt_Sstd.append(float(data[i][8]))
    network_S.append(float(data[i][9]))
    network_Sstd.append(float(data[i][11]))
    meta_S.append(float(data[i][12]))
    meta_Sstd.append(float(data[i][14]))
    meta_S2.append(float(data[i][15]))
    meta_Sstd2.append(float(data[i][16]))

timeslice= list(range(0, len(data)*3,3))

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial']})

# Create axes 
fig = plt.figure(figsize=(9,6))
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('time (min)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('S (pM/min)',fontname="Arial",fontweight="normal",fontsize="20",color='k')
#ax1.tick_params(direction='in', pad=8)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xlim([0,420])
ax1.set_ylim([0,400])
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(100)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(50)
ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
#rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
    tick.label.set_color('k')
for axis in ['bottom','left']:ax1.spines[axis].set_linewidth(2)
for axis in ['top','right']:ax1.spines[axis].set_linewidth(2)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.xaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

# Plot for heating

ax1.plot(timeslice,meal_S,linestyle='none',marker='.',markersize=4,markerfacecolor="None",markeredgecolor='k', markeredgewidth=1, alpha=0.6)
ax1.errorbar(timeslice,meal_S, yerr=meal_Sstd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.4, elinewidth=1,linewidth=1,c='k',alpha=0.6)

ax1.plot(timeslice,meta_S,linestyle='none',marker='.',markersize=4,markerfacecolor="None",markeredgecolor='red', markeredgewidth=1, alpha=0.6)
ax1.errorbar(timeslice,meta_S, yerr=meta_Sstd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.4, elinewidth=1,linewidth=1,c='red',alpha=0.6)

ax1.plot(timeslice,meta_S2,linestyle='none',marker='.',markersize=4,markerfacecolor="None",markeredgecolor='darkgreen', markeredgewidth=1, alpha=0.6)
ax1.errorbar(timeslice,meta_S2, yerr=meta_Sstd2,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.4, elinewidth=1,linewidth=1,c='darkgreen',alpha=0.6)

Sfit1 = savgol_filter(meal_S, 139, 15) # window size 51, polynomial order 3

ax1.plot(timeslice[0:50],Sfit1[0:50],linestyle='-',c='k',linewidth=1.5)
ax1.plot(timeslice[49:],meal_S[49:],linestyle='-',c='k',linewidth=1.5)

#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax1.plot((225, 250),(383,383),linestyle='-',c='k',linewidth=1.5)
ax1.plot((225, 250),(360,360),linestyle='-',c='red',linewidth=1.5)
ax1.plot((225, 250),(337,337),linestyle='-',c='green',linewidth=1.5)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.61*(left+right), 0.95*(bottom+top), 'S.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.61*(left+right), 0.89*(bottom+top), 'S.Meta without GLP-1',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.61*(left+right), 0.83*(bottom+top), 'S.Meta with GLP-1',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)

ax1.text(0.02*(left+right), 0.95*(bottom+top), 'T2D subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.88*(bottom+top), '[GLP-1] = 0.0 mM',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meta_GS_t2d_Gobs.png", transparent=True,dpi=400)
