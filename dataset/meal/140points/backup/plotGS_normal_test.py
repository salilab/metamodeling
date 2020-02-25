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

def meal_model_normal(G, DG, dt):
    Y=[None]*len(G)
    Y[0] = 0
    S=[None]*len(G)
    S[0] = 34
    for i in range(0, len(G)-1):
        #print('G',G[i])
        Y[i+1] = (1 - 0.05*dt)*Y[i] + 0.05*dt*39.6*(G[i]-5.1)
        S[i+1] = Y[i+1] + 34 + 828 * DG[i]/dt
    return S

def meal_model_t2d(G, DG, dt):
    Y=[None]*len(G)
    Y[0] = 0
    S=[None]*len(G)
    S[0] = 102.5
    for i in range(0, len(G)-1):
        #print('G',G[i])
        Y[i+1] = (1 - 0.013*dt)*Y[i] + 0.013*dt*22.5*(G[i]-9.2)
        S[i+1] = Y[i+1] + 102.5 + 445.5 * DG[i]/dt
    return S

# Read data file

data1 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_normal2.dat")
data2 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_t2d2.dat")
data3 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_normal2_DG.dat")
data4 = np.loadtxt("/Users/lipingsun/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_t2d2_DG.dat")

Gnormal=[]
Gt2d=[]
DGnormal=[]
DGt2d=[]
for i in range(0, len(data1)):
    Gnormal.append(float(data1[i][1]))
    Gt2d.append(float(data2[i][1]))
    DGnormal.append(float(data3[i][1]))
    DGt2d.append(float(data4[i][1]))
    
timeslice= list(range(0, len(data1)*3,3))

# calculation
#print(meal_model_normal(Gnormal))

S1=meal_model_normal(Gnormal,DGnormal, 3)
S2=meal_model_t2d(Gt2d,DGt2d, 3)
#sys.exit()

#----------write out the measurement----------
f1=open("data_normal_G-S.dat","w")
f2=open("data_t2d_G-S.dat","w")

for x in zip(Gnormal, S1):
    f1.write("{0}\t{1}\n".format(*x))
    
for x in zip(Gt2d, S2):
    f2.write("{0}\t{1}\n".format(*x))
    
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
ax1.set_ylabel('Gex.Meal (mM)',fontname="Arial",fontweight="normal",fontsize="20",color='chocolate')
#ax1.tick_params(direction='in', pad=8)
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
    tick.label.set_color('chocolate')
for axis in ['bottom','left']:ax1.spines[axis].set_linewidth(2)
for axis in ['top','right']:ax1.spines[axis].set_linewidth(2)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
#rcParams['xtick.direction'] = 'in'
#rcParams['ytick.direction'] = 'in'
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

#
ax2 = ax1.twinx()
ax2.yaxis.tick_right()
ax2.set_ylabel('S.Meal (pM/min)',fontname="Arial",fontweight="normal",fontsize="20",color='green')
ax2.yaxis.set_label_coords(1.1, 0.5)
ax2.set_ylim([0,400])
#ax2.tick_params(direction='in', pad=6)
xmajorLocator   = MultipleLocator(100)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
ax2.xaxis.set_major_locator(xmajorLocator)
ax2.xaxis.set_major_formatter(xmajorFormatter)
ax2.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(100)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(50)
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
ytick_locs=[0,100,200,300,400]
ytick_lbls=['0','100','200','300','400']
plt.yticks(ytick_locs,ytick_lbls)
ax2.set_yticklabels(ytick_lbls,fontname="Arial",fontweight="normal",fontsize="20",color='green')

# Plot for heating
ax1.plot(timeslice,Gnormal,linestyle='none',marker='o',markersize=4,markerfacecolor="None",markeredgecolor='chocolate', markeredgewidth=1, alpha=0.6)

ax1.plot(timeslice,Gnormal,linestyle='-',c='brown',linewidth=1.5)

ax2.plot(timeslice,S1,linestyle='none',marker='o',markersize=4,markerfacecolor="None",markeredgecolor='green', markeredgewidth=1, alpha=0.6)

Sfit1 = savgol_filter(S1, 139, 15) # window size 51, polynomial order 3
Sfit2 = savgol_filter(S2, 139, 15) # window size 51, polynomial order 3

ax2.plot(timeslice[0:50],Sfit1[0:50],linestyle='-',c='darkgreen',linewidth=1.5)
ax2.plot(timeslice[49:],S1[49:],linestyle='-',c='darkgreen',linewidth=1.5)

#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax1.plot((310, 335),(23,23),linestyle='-',c='brown',linewidth=1.5)
ax1.plot((310, 335),(21.2,21.2),linestyle='-',c='darkgreen',linewidth=1.5)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.82*(left+right), 0.95*(bottom+top), 'Gex.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.82*(left+right), 0.88*(bottom+top), 'S.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.88*(bottom+top), 'T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='r',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meal_GS_test_normal.png", transparent=True,dpi=400)
