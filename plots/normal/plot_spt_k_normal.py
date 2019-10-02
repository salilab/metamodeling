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
from math import *
import csv

# The probability density of the normal distribution
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) /np.power(2*np.pi*np.power(sig, 2.),0.5)

with open("../../bnet/meta_spt_k_normal.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()

spt_k_mu=[]
spt_k_sigma=[]
spt_k_std=[]
meta_k_mu=[]
meta_k_sigma=[]
meta_k_std=[]
glp1_k_mu=[]
glp1_k_sigma=[]
glp1_k_std=[]

for i in range(0, len(data)):
    spt_k_mu.append(float(data[i][0]))
    spt_k_sigma.append(float(data[i][1]))
    spt_k_std.append(float(data[i][2]))
    meta_k_mu.append(float(data[i][3]))
    meta_k_sigma.append(float(data[i][4]))
    meta_k_std.append(float(data[i][5]))
    glp1_k_mu.append(float(data[i][6]))
    glp1_k_sigma.append(float(data[i][7]))
    glp1_k_std.append(float(data[i][8]))
    
# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial']})

# Create axes 
fig = plt.figure(figsize=(8.5,6)) #cm
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('k.Meta',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('probability',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xlim([0,20])
ax1.set_ylim([0,1])
xmajorLocator   = MultipleLocator(6)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(3)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.2)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.1)
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
for axis in ['top','right']:ax1.spines[axis].set_linewidth(0)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'    
for line in ax1.xaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
#ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
#xtick_locs=[0.00,0.03,0.05,0.10,0.20]
#xtick_lbls=['0.00','0.03',' 0.05','0.10','0.20']
#plt.xticks(xtick_locs,xtick_lbls)
#ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
#ax1.set_yticks([0.015,0.025,0.035,0.045,0.055])
#ax1.set_axis_bgcolor('none')
#ax1.grid(True)

#prior-[10.0000]    [4.0000]    [2.0000]
#posteror-[10.2748]    [3.9993]    [1.9998]
#meta - [10.6198]    [3.9993]    [1.9998]

x=np.linspace(0, 20, 100)
for i in range(130,131):
    spt_k = gaussian(x, spt_k_mu[i], spt_k_sigma[i])
    meta_k = gaussian(x, meta_k_mu[i], meta_k_sigma[i])
    glp1_k = gaussian(x, glp1_k_mu[i], glp1_k_sigma[i])
    ax1.plot((spt_k_mu[i],spt_k_mu[i]),(0,max(spt_k)),'k',linestyle=":",linewidth=1.5)
    ax1.plot((meta_k_mu[i],meta_k_mu[i]),(0,max(meta_k)),'red',linestyle=":",linewidth=1.5)
    #ax1.plot((glp1_k_mu[i],glp1_k_mu[i]),(0,max(glp1_k)),'green',linestyle=":",linewidth=1.5)
    print(spt_k_mu[i], spt_k_std[i])
    print(meta_k_mu[i], meta_k_std[i])
    print(glp1_k_mu[i], glp1_k_std[i])

ax1.plot(x,spt_k,'k',linestyle="-",linewidth=2)
ax1.plot(x,meta_k,'red',linestyle="-",linewidth=2)
#ax1.plot(x,glp1_k,'green',linestyle="-",linewidth=2)

#ax1.plot((15, 17.3),(0.96,0.96),linestyle='-',c='k',linewidth=1.5)
#ax1.plot((15, 17.3),(0.9,0.9),linestyle='-',c='red',linewidth=1.5)
#ax1.plot((15, 17.3),(0.84,0.84),linestyle='-',c='green',linewidth=1.5)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
#ax1.text(0.6*(left+right), 0.95*(bottom+top), 'k.SPT',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.6*(left+right), 0.89*(bottom+top), 'k.Meta without GLP-1',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.6*(left+right), 0.83*(bottom+top), 'k.Meta with GLP-1',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.84*(left+right), 0.88*(bottom+top), '$V_{ISG}$ - Meta',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meta_spt_k_normal.png",dpi=300)






