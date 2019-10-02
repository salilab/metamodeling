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

data = np.loadtxt("meal_exp1_normal2_Gin_dt1_sigmoid.dat")
DG = data[:,1]
G = data[:,2]

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
ax1.set_ylabel('Gintake.Meal (mM)',fontname="Arial",fontweight="normal",fontsize="20",color='k')
#ax1.tick_params(direction='in', pad=8)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.08, 0.5)
ax1.set_xlim([0,420])
ax1.set_ylim([0,12])
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(3)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(1.5)
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
#ax1.errorbar(timeslice,G, yerr=Gstd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.4, elinewidth=1,linewidth=1,c='chocolate',alpha=0.6)
ax1.plot(timeslice,G,linestyle='-',c='k',linewidth=1.5)
ax1.plot(timeslice,DG,linestyle='-',c='r',linewidth=1.5)
#ax2.plot(timeslice,I,linestyle='-',marker='.',markersize=4,markerfacecolor="None",markeredgecolor='green', markeredgewidth=1, alpha=0.6)
#ax2.errorbar(timeslice,I, yerr=Istd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.4, elinewidth=1,linewidth=1,c='green',alpha=0.6)

#Ifit1 = savgol_filter(I[0:], 419, 20) # window size 51, polynomial order 3
#ax2.plot(timeslice[0:],Ifit1,linestyle='-',c='darkgreen',linewidth=1.5)
#ax2.plot(timeslice[49:],I[49:],linestyle='-',c='darkgreen',linewidth=1.5)

#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
#ax1.plot((310, 335),(23,23),linestyle='-',c='brown',linewidth=1.5)
#ax1.plot((310, 335),(21.2,21.2),linestyle='-',c='darkgreen',linewidth=1.5)
#ax1.plot((60,60),(0, 30),linestyle=':',c='k',linewidth=1.5)

ax1.plot((60, 60),(0, 30),linestyle='--',c='k',linewidth=1)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
#ax1.text(0.82*(left+right), 0.95*(bottom+top), 'Gex.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.82*(left+right), 0.88*(bottom+top), 'I.Meal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.88*(bottom+top), 'Evidence = Gex.obs, DGex.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("data_normal2_Gin_dt1_sigmoid.png", transparent=True,dpi=400)
