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

# define function

def meal_model_normal(G):
    Y=[None]*len(G)
    Y[0] = 0
    Spo=[None]*len(G)
    Spo[0] = 34
    for i in range(0, len(G)-1):
        #print('G',G[i])
        Y[i+1] = (1 - 0.05*1)*Y[i] + 0.05*1*39.6*(G[i]-5.1)
        if G[i+1]-G[i] > 0:
            #print(G[i+1]-G[i])
            Spo[i+1] = (Y[i+1] + 34 + 828 * ((G[i+1]-G[i])/1))
        else:
            Spo[i+1] = (Y[i+1] + 34)
    return Spo

def meal_model_t2b(G):
    Y=[None]*len(G)
    Y[0] = 0
    Spo=[None]*len(G)
    Spo[0] = 102.5
    for i in range(0, len(G)-1):
        #print('G',G[i])
        Y[i+1] = (1 - 0.013*1)*Y[i] + 0.013*1*22.5*(G[i]-9.2)
        if G[i+1]-G[i] > 0:
            Spo[i+1] = (Y[i+1] + 102.5 + 445.5 * ((G[i+1]-G[i])/1))
        else:
            Spo[i+1] = (Y[i+1] + 102.5)
    return Spo

# Read data file

data1 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modeling2/dataset/meal/420points/meal_exp1_normal2.dat")
data2 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modeling2/dataset/meal/420points/meal_exp1_t2b2.dat")
   
Gnormal=[]
Gt2b=[]
for i in range(0, len(data1)):
    Gnormal.append(float(data1[i][1]))
    Gt2b.append(float(data2[i][1]))

timeslice= list(range(0, len(data1)))

# calculation
#print(meal_model_normal(Gnormal))

I1=meal_model_normal(Gnormal)
I2=meal_model_t2b(Gt2b)
#sys.exit()
print(I1)
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
ax1.set_xlim([0,420])
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
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax2.xaxis.set_major_locator(xmajorLocator)
ax2.xaxis.set_major_formatter(xmajorFormatter)
ax2.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(80)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(40)
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

ax1.plot(timeslice,Gnormal,marker='.',linestyle='none',c='k',markersize=2)
ax1.plot(timeslice,Gt2b,marker='.',linestyle='none',c='r',markersize=2)
Gfit1 = savgol_filter(Gnormal, 419, 10) # window size 51, polynomial order 3
Gfit2 = savgol_filter(Gt2b, 419, 10) # window size 51, polynomial order 3
#ax1.plot(timeslice,Gfit1,linestyle='-',c='k',linewidth=2)
#ax1.plot(timeslice,Gfit2,linestyle='-',c='r',linewidth=2)

Ifit1 = savgol_filter(I1, 419, 12) # window size 51, polynomial order 3
Ifit2 = savgol_filter(I2, 419, 12) # window size 51, polynomial order 3

ax2.plot(timeslice,I1,marker='o',linestyle='none',c='green',markersize=2)
ax2.plot(timeslice,I2,marker='o',linestyle='none',c='blue',markersize=2)

ax2.plot(timeslice[200:],I1[200:],linestyle=':',c='k',linewidth=2)
ax2.plot(timeslice[200:],I2[200:],linestyle=':',c='r',linewidth=2)
ax2.plot(timeslice[0:200],Ifit1[0:200],linestyle=':',c='k',linewidth=2)
ax2.plot(timeslice[0:200],Ifit2[0:200],linestyle=':',c='r',linewidth=2)
#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax1.plot((340, 365),(23,23),linestyle='-',c='k',linewidth=2)
ax1.plot((340, 365),(21.2,21.2),linestyle='--',c='k',linewidth=2)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.88*(left+right), 0.95*(bottom+top), 'G.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.88*(left+right), 0.88*(bottom+top), 'I.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='r',transform=ax1.transAxes)
ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.02*(left+right), 0.88*(bottom+top), 'T2B',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='r',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meal_GI_test.jpg",dpi=400)
