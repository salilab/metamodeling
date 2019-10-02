#!/usr/bin/python

# Import module for plot "matplotlib"%
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
with open("../../bnet/meta_spt_normal_variables.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()

with open("../../bnet/meta_spt_t2d_variables.txt", 'rb') as kfile:
    next(kfile)
    data2 = list(csv.reader(kfile))
    kfile.close()
    
k_prior_mu=[]
k_prior_sigma=[]
k_posterior_mu=[]
k_posterior_sigma=[]
for i in range(0, len(data)):
    k_prior_mu.append(float(data[i][5]))
    k_prior_sigma.append(float(data[i][6]))
    k_posterior_mu.append(float(data[i][26]))
    k_posterior_sigma.append(float(data[i][27]))
        
k_prior_mu2=[]
k_prior_sigma2=[]
k_posterior_mu2=[]
k_posterior_sigma2=[]
for i in range(0, len(data)):
    k_prior_mu2.append(float(data2[i][3]))
    k_prior_sigma2.append(float(data2[i][4]))
    k_posterior_mu2.append(float(data2[i][24]))
    k_posterior_sigma2.append(float(data2[i][25]))
    
print(k_posterior_mu)    
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
ax1.set_xlabel('Npatch.SPT',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('Probability density function',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_xlim([0,20])
ax1.set_ylim([0,0.2])
xmajorLocator   = MultipleLocator(4)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(2)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.05)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.025)
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

print(k_prior_mu[99], k_prior_sigma[99])
print(k_posterior_mu[99], k_posterior_sigma[99])

k_prior = gaussian(x, k_prior_mu[99], k_prior_sigma[99])
k_posterior = gaussian(x, k_posterior_mu[99], k_posterior_sigma[99])

k_prior2 = gaussian(x, k_prior_mu2[99], k_prior_sigma2[99])
k_posterior2 = gaussian(x, k_posterior_mu2[99], k_posterior_sigma2[99])

ax1.plot(x,k_prior,'k',linestyle="-",linewidth=2)
ax1.plot(x,k_posterior,'red',linestyle="-",linewidth=2)

ax1.plot((k_prior_mu[99],k_prior_mu[99]),(0,max(k_prior)),'k',linestyle=":",linewidth=1)
ax1.plot((k_posterior_mu[99],k_posterior_mu[99]),(0,max(k_posterior)),'red',linestyle=":",linewidth=1)

#ax1.plot(x,k_prior2,'green',linestyle="-",linewidth=2)
ax1.plot(x,k_posterior2,'darkorange',linestyle="-",linewidth=2)

#ax1.plot((k_prior_mu2[99],k_prior_mu2[99]),(0,max(k_prior2)),'green',linestyle=":",linewidth=1)
ax1.plot((k_posterior_mu2[99],k_posterior_mu2[99]),(0,max(k_posterior2)),'darkorange',linestyle=":",linewidth=1)

ax1.plot((13.3, 14.5),(0.193, 0.193),linestyle='-',c='k',linewidth=1.5)
ax1.plot((13.3, 14.5),(0.182, 0.182),linestyle='-',c='red',linewidth=1.5)
ax1.plot((13.3, 14.5),(0.171, 0.171),linestyle='-',c='darkorange',linewidth=1.5)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.75*(left+right), 0.96*(bottom+top), 'Npatch.SPT - Normal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.75*(left+right), 0.91*(bottom+top), 'Npatch.Meta - Normal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.75*(left+right), 0.86*(bottom+top), 'Npatch.Meta - T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.84*(left+right), 0.88*(bottom+top), '$V_{ISG}$ - Meta',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meta_spt_Npatch_normal.png",dpi=300)






