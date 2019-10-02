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

# Read data file
with open("../normal/k_spt_prior-posterior.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()
    
k_prior_mu=[]
k_prior_sigma=[]
k_posterior_mu=[]
k_posterior_sigma=[]
for i in range(0, len(data)):
    k_prior_mu.append(float(data[i][0]))
    k_prior_sigma.append(float(data[i][1]))
    k_posterior_mu.append(float(data[i][3]))
    k_posterior_sigma.append(float(data[i][4]))
    
#k2 = np.loadtxt("output/random_pbc_NE-RDF1-2.xvg")

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial']})

# Create axes 
fig = plt.figure(figsize=(8.5,6)) #cm
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('k.SPT',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('probability',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_xlim([0,20])
ax1.set_ylim([0,0.4])
xmajorLocator   = MultipleLocator(4)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(2)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.1)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.05)
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
#ax1.plot((0.2,0.2),(0,1),'grey',linestyle=":",linewidth=1.5)
#10.2091+-0.6346

x=np.linspace(0, 20, 400)

print(k_prior_mu[-1], k_prior_sigma[-1])
print(k_posterior_mu[-1], k_posterior_sigma[-1])

k_prior = gaussian(x, k_prior_mu[-1], k_prior_sigma[-1])
k_posterior = gaussian(x, k_posterior_mu[-1], k_posterior_sigma[-1])

ax1.plot(x,k_prior,'k',linestyle="-",linewidth=2)
ax1.plot(x,k_posterior,'red',linestyle="-",linewidth=2)

ax1.plot((k_prior_mu[-1],k_prior_mu[-1]),(0,max(k_prior)),'k',linestyle=":",linewidth=2)
ax1.plot((k_posterior_mu[-1],k_posterior_mu[-1]),(0,max(k_posterior)),'red',linestyle=":",linewidth=2)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.96*(left+right), 0.96*(bottom+top), 'with no GLP-1R: 10.0 $\pm$ 2.0',horizontalalignment='right',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.96*(left+right), 0.88*(bottom+top), 'with GLP-1R: 9.5 $\pm$ 1.7',horizontalalignment='right',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
#ax1.text(0.84*(left+right), 0.88*(bottom+top), '$V_{ISG}$ - Meta',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("normal_k_meta_posterior3-4.png",dpi=300)






