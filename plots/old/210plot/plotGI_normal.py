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

# Read data file

with open("../../meal_G-Spo_normal.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()
    
G=[]
Gstd=[]
Spo=[]
Spostd=[]
I=[]
for i in range(0, len(data)):
    G.append(float(data[i][0]))
    Gstd.append(float(data[i][1]))
    Spo.append(float(data[i][2]))
    Spostd.append(float(data[i][3]))
    I.append(float(data[i][4]))
    
timeslice= list(range(0, len(data)))

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
ax1.set_ylabel('G.Meal (mM)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=8)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.08, 0.5)
ax1.set_xlim([0,84])
ax1.set_ylim([0,24])
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(6)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(3)
ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
for axis in ['bottom','left']:ax1.spines[axis].set_linewidth(2)
for axis in ['top','right']:ax1.spines[axis].set_linewidth(2)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

#
ax2 = ax1.twinx()
ax2.yaxis.tick_right()
ax2.set_ylabel('Spo.Meal (pM/min)',fontname="Arial",fontweight="normal",fontsize="20",color='red')
ax2.yaxis.set_label_coords(1.1, 0.5)
ax2.set_ylim([0,400])
ax2.tick_params(direction='in', pad=6)
xmajorLocator   = MultipleLocator(50)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(25)
ax2.xaxis.set_major_locator(xmajorLocator)
ax2.xaxis.set_major_formatter(xmajorFormatter)
ax2.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(100)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(50)
ax2.yaxis.set_major_locator(ymajorLocator)
ax2.yaxis.set_major_formatter(ymajorFormatter)
ax2.yaxis.set_minor_locator(yminorLocator)
rcParams['ytick.direction'] = 'in'
for line in ax2.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax2.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
ytick_locs=[0,100,200,300,400]
ytick_lbls=['0','100','200','300','400']
plt.yticks(ytick_locs,ytick_lbls)
ax2.set_yticklabels(ytick_lbls,fontname="Arial",fontweight="normal",fontsize="20",color='r')

# Plot for heating
ax1.plot(timeslice,G,linestyle='-',c='k',linewidth=2)
ax1.errorbar(timeslice,G, yerr=Gstd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=0.2, elinewidth=0.2,linewidth=2,c='k')

ax2.plot(timeslice,Spo,linestyle='-',c='r',linewidth=2)
ax2.errorbar(timeslice,Spo, yerr=Spostd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=0.2, elinewidth=0.2,linewidth=2,c='r')

#ax2.plot(timeslice,I,linestyle='-',c='r',linewidth=2)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.8*(left+right), 0.95*(bottom+top), 'G.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.8*(left+right), 0.88*(bottom+top), 'Spo.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='r',transform=ax1.transAxes)
ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meal_GI_normal.jpg",dpi=400)
