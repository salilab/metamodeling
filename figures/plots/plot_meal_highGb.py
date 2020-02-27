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
with open("../../bnet/meta_meal_normal_variables_highGb.txt", 'rb') as kfile:
    next(kfile)
    data = list(csv.reader(kfile))
    kfile.close()

with open("../../bnet/meta_meal_t2d_variables.txt", 'rb') as kfile:
    next(kfile)
    data2 = list(csv.reader(kfile))
    kfile.close()
    
Gb_prior_mu=[]
Gb_prior_sigma=[]
Gb_posterior_mu=[]
Gb_posterior_sigma=[]
for i in range(0, len(data)):
    Gb_prior_mu.append(float(data[i][6]))
    Gb_prior_sigma.append(float(data[i][7]))
    Gb_posterior_mu.append(float(data[i][18]))
    Gb_posterior_sigma.append(float(data[i][19]))
        
Gb_prior_mu2=[]
Gb_prior_sigma2=[]
Gb_posterior_mu2=[]
Gb_posterior_sigma2=[]
for i in range(0, len(data)):
    Gb_prior_mu2.append(float(data2[i][6]))
    Gb_prior_sigma2.append(float(data2[i][7]))
    Gb_posterior_mu2.append(float(data2[i][18]))
    Gb_posterior_sigma2.append(float(data2[i][19]))
    
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
ax1.set_xlabel('G$_{basal}$ (mM)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('Probability density function',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.12, 0.5)
ax1.set_xlim([4,10])
ax1.set_ylim([0,8])
xmajorLocator   = MultipleLocator(1)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(0.5)
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
ax1.plot((5.111,5.111),(0,10),'k',linestyle="--",linewidth=1.5)
#10.2091+-0.6346

x=np.linspace(0, 20, 1000)

print(Gb_prior_mu[400], Gb_prior_sigma[400])
print(Gb_posterior_mu[400], Gb_posterior_sigma[400])

Gb_prior = gaussian(x, Gb_prior_mu[400], Gb_prior_sigma[400])
Gb_posterior = gaussian(x, Gb_posterior_mu[400], Gb_posterior_sigma[400])

Gb_prior2 = gaussian(x, Gb_prior_mu2[400], Gb_prior_sigma2[400])
Gb_posterior2 = gaussian(x, Gb_posterior_mu2[400], Gb_posterior_sigma2[400])

ax1.plot(x,Gb_prior,'k',linestyle="-",linewidth=2)
#ax1.plot(x,Gb_posterior,'blue',linestyle="-",linewidth=2)

ax1.plot((Gb_prior_mu[400],Gb_prior_mu[400]),(0,max(Gb_prior)),'k',linestyle=":",linewidth=1)
#ax1.plot((Gb_posterior_mu[400],Gb_posterior_mu[400]),(0,max(Gb_posterior)),'blue',linestyle=":",linewidth=1)

#ax1.plot(x,Gb_prior2,'brown',linestyle="-",linewidth=2)
#ax1.plot(x,Gb_posterior2,'darkorange',linestyle="-",linewidth=2)

#ax1.plot((Gb_prior_mu2[400],Gb_prior_mu2[400]),(0,max(Gb_prior2)),'brown',linestyle=":",linewidth=1)
#ax1.plot((Gb_posterior_mu2[400],Gb_posterior_mu2[400]),(0,max(Gb_posterior2)),'darkorange',linestyle=":",linewidth=1)

#ax1.plot((7.8, 8.4),(2.4, 2.4),linestyle='-',c='k',linewidth=1.5)
#ax1.plot((7.8, 8.4),(2.2, 2.2),linestyle='-',c='brown',linewidth=1.5)

p=patches.FancyArrowPatch((4.6,0.75),(3,0.75),arrowstyle='->',linewidth=2,connectionstyle='arc3',mutation_scale=20,alpha=1)
p.set_edgecolor('k')
p.set_facecolor('k')
p.set_zorder(10)
#ax1.add_patch(p)

p=patches.FancyArrowPatch((8.8,0.75),(7,0.75),arrowstyle='->',linewidth=2,connectionstyle='arc3',mutation_scale=20,alpha=1)
p.set_edgecolor('r')
p.set_facecolor('r')
p.set_zorder(10)
#ax1.add_patch(p)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.21*(left+right), 0.95*(bottom+top), 'Literature \n value',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.71*(left+right), 0.88*(bottom+top), 'Gb.Meal - T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.32*(left+right), 0.91*(bottom+top), 'G$_b$.Meta - Normal',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.75*(left+right), 0.91*(bottom+top), 'G$_b$.Meta - T2D',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.84*(left+right), 0.88*(bottom+top), '$V_{ISG}$ - Meta',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meal_Gb_highGb.png",dpi=300)






