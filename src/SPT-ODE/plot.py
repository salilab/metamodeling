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

with open("k.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()
    
sptk=[]
sptkstd=[]
metak=[]
metakstd=[]
for i in range(0, len(data)):
    sptk.append(float(data[i][0]))
    sptkstd.append(float(data[i][1]))
    metak.append(float(data[i][2]))
    metakstd.append(float(data[i][3]))
del data


with open("h.txt", 'rb') as hfile:
    next(hfile)
    data = list(csv.reader(hfile))
    kfile.close()
    
odeh=[]
odehstd=[]
metah=[]
metahstd=[]
for i in range(0, len(data)):
    odeh.append(float(data[i][0]))
    odehstd.append(float(data[i][1]))
    metah.append(float(data[i][2]))
    metahstd.append(float(data[i][3]))
del data

# plasma insulin concentration after a meal - pmol/l
insulin=[]
for i in range(0,50):
    j = i*10.0 
    insulin.append(j)

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial']})

# Create axes 
fig = plt.figure(figsize=(9,6))
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('Plasma insulin',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('SPT.k',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=8)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xlim([0,500])
ax1.set_ylim([8,16])
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(2)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(1)
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
ax2.set_ylabel('ODE.h',fontname="Arial",fontweight="normal",fontsize="20")
#ax2.yaxis.set_label_coords(-0.1, 0.5)
ax2.set_ylim([5.5,6.5])
ax2.tick_params(direction='in', pad=6)
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax2.xaxis.set_major_locator(xmajorLocator)
ax2.xaxis.set_major_formatter(xmajorFormatter)
ax2.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.25)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.125)
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
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%2f'))
ytick_locs=[5.50,5.75,6.00,6.25,6.50]
ytick_lbls=['5.50','5.75','6.00','6.25','6.50']
plt.yticks(ytick_locs,ytick_lbls)
ax2.set_yticklabels(ytick_lbls,fontname="Arial",fontweight="normal",fontsize="20",color='r')

#ax1.set_axis_bgcolor('none')
#ax1.grid(True)
#ax1.plot((0,0),(0,350),'grey',linestyle=":",linewidth=1.5)


# Plot for heating
ax1.plot(insulin,sptk,linestyle='--',c='k',linewidth=2)
#ax1.errorbar(insulin,sptk, yerr=sptkstd,marker='o',linestyle='none',markersize=6, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='k')
ax1.plot(insulin,metak,linestyle='-',c='k',linewidth=2)
#ax1.errorbar(insulin,metak, yerr=metakstd,marker='o',linestyle='none',markersize=6, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='blue')

ax2.plot(insulin,odeh,linestyle='--',c='r',linewidth=2)
#ax2.errorbar(insulin,odeh, yerr=odehstd,marker='o',linestyle='none',markersize=6, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax2.plot(insulin,metah,linestyle='-',c='r',linewidth=2)

#ax2.errorbar(insulin,metah, yerr=metahstd,marker='o',linestyle='none',markersize=6, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='orange')
ax1.plot((10,40),(15.6,15.6),'k',linestyle="-",linewidth=2)
ax1.plot((10,40),(15,15),'k',linestyle="--",linewidth=2)

ax1.plot((295,325),(15.6,15.6),'r',linestyle="-",linewidth=2)
ax1.plot((295,325),(15,15),'r',linestyle="--",linewidth=2)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.09*(left+right), 0.95*(bottom+top), 'SPT.k posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.09*(left+right), 0.88*(bottom+top), 'meta.k posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.66*(left+right), 0.95*(bottom+top), 'ODE.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("h-k_insulin.jpg",dpi=400)