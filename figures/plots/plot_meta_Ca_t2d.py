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
with open("../../bnet/meta_network_normal_variables.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()

with open("../../bnet/meta_network_t2d_variables.txt", 'rb') as kfile:
    next(kfile)
    data2 = list(csv.reader(kfile))
    kfile.close()
    
k_prior_mu=[]
k_prior_sigma=[]
k_posterior_mu=[]
k_posterior_sigma=[]
for i in range(0, len(data)):
    k_prior_mu.append(float(data[i][6])*10000)
    k_prior_sigma.append(float(data[i][7])*10000.0)
    k_posterior_mu.append(float(data[i][21])*10000)
    k_posterior_sigma.append(float(data[i][23])*10000.0)
        
k_prior_mu2=[]
k_prior_sigma2=[]
k_posterior_mu2=[]
k_posterior_sigma2=[]
for i in range(0, len(data)):
    k_prior_mu2.append(float(data2[i][6])*10000)
    k_prior_sigma2.append(float(data2[i][8])*10000.0)
    k_posterior_mu2.append(float(data2[i][21])*10000)
    k_posterior_sigma2.append(float(data2[i][23])*10000.0)
    
#print(k_posterior_mu)    
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
ax1.set_xlabel('Ca$^{2+}$ ( 10$^{-4}$ mM)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('Probability density function',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_xlim([0.6,1.4])
ax1.set_ylim([0,15])
xmajorLocator   = MultipleLocator(0.2)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.1)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(3)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(1.5)
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

x=np.linspace(0, 2, 2000)

print(k_posterior_mu[499], k_posterior_sigma[499])
print(k_posterior_mu2[499], k_posterior_sigma2[499])

k_prior = gaussian(x, k_prior_mu[499], k_prior_sigma[499])
k_posterior = gaussian(x, k_posterior_mu[499], k_posterior_sigma[499])

k_prior2 = gaussian(x, k_prior_mu2[499], k_prior_sigma2[499])
k_posterior2 = gaussian(x, k_posterior_mu2[499], k_posterior_sigma2[499])

#ax1.plot(x,k_prior,'k',linestyle="-",linewidth=2)
#ax1.plot(x,k_posterior,'k',linestyle="-",linewidth=2)

#ax1.plot((k_prior_mu[499],k_prior_mu[499]),(0,max(k_prior)),'k',linestyle=":",linewidth=1)
#ax1.plot((k_posterior_mu[499],k_posterior_mu[499]),(0,max(k_posterior)),'k',linestyle=":",linewidth=1)

#ax1.plot(x,k_prior2,'green',linestyle="-",linewidth=2)
ax1.plot(x,k_posterior2,'r',linestyle="-",linewidth=2)

#ax1.plot((k_prior_mu2[499],k_prior_mu2[499]),(0,max(k_prior2)),'green',linestyle=":",linewidth=1)
ax1.plot((k_posterior_mu2[499],k_posterior_mu2[499]),(0,max(k_posterior2)),'r',linestyle=":",linewidth=1)

#ax1.plot((15, 16),(0.192, 0.192),linestyle='-',c='k',linewidth=1.5)
#ax1.plot((15, 16),(0.178, 0.178),linestyle='-',c='red',linewidth=1.5)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.02*(left+right), 0.95*(bottom+top), 'T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.82*(left+right), 0.95*(bottom+top), 'Ca.SPT',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.82*(left+right), 0.88*(bottom+top), 'k.Meta',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.75*(left+right), 0.86*(bottom+top), 'k.Meta - T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.84*(left+right), 0.88*(bottom+top), '$V_{ISG}$ - Meta',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meta_netework_Ca_t2d.png",dpi=300)





