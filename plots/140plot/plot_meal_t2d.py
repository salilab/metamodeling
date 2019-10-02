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
with open("../../bnet/meal_t2d.txt", 'rb') as kfile:
    next(kfile)
    meal = list(csv.reader(kfile))
    kfile.close()

meal_G=[]
meal_Gstd=[]
meal_S=[]
meal_Sstd=[]
meal_I=[]
meal_Istd=[]

for i in range(0, len(meal)):
    meal_G.append(float(meal[i][0]))
    meal_Gstd.append(float(meal[i][2]))
    meal_I.append(float(meal[i][6]))
    meal_Istd.append(float(meal[i][8]))
    
timeslice= list(range(0, len(meal)))

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
ax1.set_ylabel('G.Meal (mM)',fontname="Arial",fontweight="normal",fontsize="20",color='chocolate')
#ax1.tick_params(direction='in', pad=8)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.08, 0.5)
ax1.set_xlim([0,420])
ax1.set_ylim([0,30])
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(10)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(5)
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
    tick.label.set_color('chocolate')
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

#
ax2 = ax1.twinx()
ax2.yaxis.tick_right()
ax2.set_ylabel('I.Meal (pM)',fontname="Arial",fontweight="normal",fontsize="20",color='darkgreen')
ax2.yaxis.set_label_coords(1.1, 0.5)
ax2.set_ylim([0,400])
#ax2.tick_params(direction='in', pad=6)
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax2.xaxis.set_major_locator(xmajorLocator)
ax2.xaxis.set_major_formatter(xmajorFormatter)
ax2.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(200)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(100)
ax2.yaxis.set_major_locator(ymajorLocator)
ax2.yaxis.set_major_formatter(ymajorFormatter)
ax2.yaxis.set_minor_locator(yminorLocator)
#rcParams['ytick.direction'] = 'in'
for line in ax2.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax2.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax2.yaxis.set_ticks_position('right')
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
ytick_locs=[0,200,400,600]
ytick_lbls=['0','200','400','600']
plt.yticks(ytick_locs,ytick_lbls)
ax2.set_yticklabels(ytick_lbls,fontname="Arial",fontweight="normal",fontsize="20",color='darkgreen')

# Plot for heating
ax1.errorbar(timeslice,meal_G, yerr=meal_Gstd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='chocolate',alpha=0.6)
ax1.plot(timeslice,meal_G,linestyle='-',c='brown',linewidth=1)

#ax2.plot(timeslice,I,linestyle='-',marker='.',markersize=4,markerfacecolor="None",markeredgecolor='green', markeredgewidth=1, alpha=0.6)
ax2.errorbar(timeslice,meal_I, yerr=meal_Istd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='green',alpha=0.6)
ax2.plot(timeslice,meal_I,linestyle='-',c='darkgreen',linewidth=1, alpha=0.6)

#Ifit1 = savgol_filter(I[0:], 419, 20) # window size 51, polynomial order 3
#ax2.plot(timeslice[0:],Ifit1,linestyle='-',c='darkgreen',linewidth=1.5)
#ax2.plot(timeslice[49:],I[49:],linestyle='-',c='darkgreen',linewidth=1.5)

#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax1.plot((327, 340),(28.5,28.5),linestyle='-',c='brown',linewidth=1.5)
ax1.plot((327, 340),(26.5,26.5),linestyle='-',c='darkgreen',linewidth=1.5)
#ax1.plot((60,60),(0, 30),linestyle=':',c='k',linewidth=1.5)

ax1.plot((0, 420),(meal_G[0],meal_G[0]),linestyle=':',c='brown',linewidth=1)
ax2.plot((0, 420),(meal_I[0], meal_I[0]),linestyle=':',c='darkgreen',linewidth=1)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.02*(left+right), 0.95*(bottom+top), 'T2D subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.82*(left+right), 0.95*(bottom+top), 'G.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.82*(left+right), 0.88*(bottom+top), 'I.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.88*(bottom+top), 'Evidence = Gex.obs, DGex.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meal_t2d.png", transparent=True,dpi=400)
