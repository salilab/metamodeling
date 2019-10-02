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
with open("../../bnet/meta_meal_t2d.txt", 'rb') as kfile:
    next(kfile)
    meal = list(csv.reader(kfile))
    kfile.close()
    
with open("../../bnet/meta_meal_t2d_glp1_medium.txt", 'rb') as kfile:
    next(kfile)
    meal1 = list(csv.reader(kfile))
    kfile.close()
    
with open("../../bnet/meta_meal_t2d_glp1_high.txt", 'rb') as kfile:
    next(kfile)
    meal2 = list(csv.reader(kfile))
    kfile.close()

meal_G=[]
meal_Gstd=[]
meal_I=[]
meal_Istd=[]
meta_meal_G=[]
meta_meal_Gstd=[]
meta_meal_I=[]
meta_meal_Istd=[]

meal_G_medium=[]
meal_Gstd_medium=[]
meal_I_medium=[]
meal_Istd_medium=[]
meta_meal_G_medium=[]
meta_meal_Gstd_medium=[]
meta_meal_I_medium=[]
meta_meal_Istd_medium=[]

meal_G_high=[]
meal_Gstd_high=[]
meal_I_high=[]
meal_Istd_high=[]
meta_meal_G_high=[]
meta_meal_Gstd_high=[]
meta_meal_I_high=[]
meta_meal_Istd_high=[]

for i in range(0, len(meal1)):
    meal_G.append(float(meal[i][0]))
    meal_Gstd.append(float(meal[i][2]))
    meal_I.append(float(meal[i][9]))
    meal_Istd.append(float(meal[i][11]))
    meta_meal_G.append(float(meal[i][12]))
    meta_meal_Gstd.append(float(meal[i][14]))
    meta_meal_I.append(float(meal[i][21]))
    meta_meal_Istd.append(float(meal[i][23]))
    
    meal_G_medium.append(float(meal1[i][0]))
    meal_Gstd_medium.append(float(meal1[i][2]))
    meal_I_medium.append(float(meal1[i][9]))
    meal_Istd_medium.append(float(meal1[i][11]))
    meta_meal_G_medium.append(float(meal1[i][12]))
    meta_meal_Gstd_medium.append(float(meal1[i][14]))
    meta_meal_I_medium.append(float(meal1[i][21]))
    meta_meal_Istd_medium.append(float(meal1[i][23]))
    
    meal_G_high.append(float(meal2[i][0]))
    meal_Gstd_high.append(float(meal2[i][2]))
    meal_I_high.append(float(meal2[i][9]))
    meal_Istd_high.append(float(meal2[i][11]))
    meta_meal_G_high.append(float(meal2[i][12]))
    meta_meal_Gstd_high.append(float(meal2[i][14]))
    meta_meal_I_high.append(float(meal2[i][21]))
    meta_meal_Istd_high.append(float(meal2[i][23]))

timeslice= list(range(0, len(meal1)))

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
ax1.set_ylabel('G (mM)',fontname="Arial",fontweight="normal",fontsize="20",color='chocolate')
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
ax2.set_ylabel('I (pM)',fontname="Arial",fontweight="normal",fontsize="20",color='darkgreen')
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
#ax1.errorbar(timeslice,meal_G_medium, yerr=meal_Gstd_medium,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='grey',alpha=0.6)
#ax1.plot(timeslice,meal_G_medium,linestyle='-',c='grey',linewidth=1)

#ax2.errorbar(timeslice,meal_I_medium, yerr=meal_Istd_medium,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='k',alpha=0.6)
#ax2.plot(timeslice,meal_I_medium,linestyle='-',c='k',linewidth=1, alpha=0.6)

ax1.plot((0, 420),(9.167, 9.167),linestyle='-',c='chocolate',linewidth=1)
ax2.plot((0, 420),(52, 52),linestyle='-',c='darkgreen',linewidth=1)

ax1.errorbar(timeslice,meta_meal_G, yerr=meta_meal_Gstd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='chocolate',alpha=0.6)
ax1.plot(timeslice,meta_meal_G,linestyle='-',c='k',linewidth=1)

ax2.errorbar(timeslice,meta_meal_I, yerr=meta_meal_Istd,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='darkgreen',alpha=0.6)
ax2.plot(timeslice,meta_meal_I,linestyle='-',c='k',linewidth=1, alpha=0.6)

ax1.errorbar(timeslice,meta_meal_G_medium, yerr=meta_meal_Gstd_medium,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='darkorange',alpha=0.6)
ax1.plot(timeslice,meta_meal_G_medium,linestyle='-',c='darkorange',linewidth=1)

ax2.errorbar(timeslice,meta_meal_I_medium, yerr=meta_meal_Istd_medium,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='lime',alpha=0.6)
ax2.plot(timeslice,meta_meal_I_medium,linestyle='-',c='lime',linewidth=1, alpha=0.6)

ax1.errorbar(timeslice,meta_meal_G_high, yerr=meta_meal_Gstd_high,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='gold',alpha=0.6)
ax1.plot(timeslice,meta_meal_G_high,linestyle='-',c='gold',linewidth=1)

ax2.errorbar(timeslice,meta_meal_I_high, yerr=meta_meal_Istd_high,marker='.',linestyle='none',markersize=2, fmt='-',capsize=1.2, elinewidth=0.5,linewidth=0.5,c='greenyellow',alpha=0.6)
ax2.plot(timeslice,meta_meal_I_high,linestyle='-',c='greenyellow',linewidth=1, alpha=0.6)

#Ifit1 = savgol_filter(I[0:], 419, 20) # window size 51, polynomial order 3
#ax2.plot(timeslice[0:],Ifit1,linestyle='-',c='darkgreen',linewidth=1.5)
#ax2.plot(timeslice[49:],I[49:],linestyle='-',c='darkgreen',linewidth=1.5)

#ax1.errorbar(timeslice,I, yerr=Istd,marker='o',linestyle='none',markersize=2, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='r')
ax1.plot((129, 142),(29.1,29.1),linestyle='-',c='chocolate',linewidth=1.5)
ax1.plot((129, 142),(27.5,27.5),linestyle='-',c='darkorange',linewidth=1.5)
ax1.plot((129, 142),(25.9,25.9),linestyle='-',c='gold',linewidth=1.5)
ax1.plot((280, 293),(29.1,29.1),linestyle='-',c='darkgreen',linewidth=1.5)
ax1.plot((280, 293),(27.5,27.5),linestyle='-',c='lime',linewidth=1.5)
ax1.plot((280, 293),(25.9,25.9),linestyle='-',c='greenyellow',linewidth=1.5)


p=patches.FancyArrowPatch((55,355),(55,515),arrowstyle='->',linewidth=2,connectionstyle='arc3',mutation_scale=20,alpha=1)
p.set_edgecolor('k')
p.set_facecolor('k')
p.set_zorder(10)
ax2.add_patch(p)

p=patches.FancyArrowPatch((250,230),(250,120),arrowstyle='->',linewidth=2,connectionstyle='arc3',mutation_scale=20,alpha=1)
p.set_edgecolor('k')
p.set_facecolor('k')
p.set_zorder(10)
ax2.add_patch(p)
# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.35*(left+right), 0.96*(bottom+top), 'G.Meta',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.71*(left+right), 0.96*(bottom+top), 'I.Meta',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.35*(left+right), 0.91*(bottom+top), 'G.Meta - Glp1+',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.71*(left+right), 0.91*(bottom+top), 'I.Meta - Glp1+',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.35*(left+right), 0.85*(bottom+top), 'G.Meta - Glp1++',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
ax1.text(0.71*(left+right), 0.85*(bottom+top), 'I.Meta - Glp1++',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.95*(bottom+top), 'Normal subject',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="bold",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.88*(bottom+top), 'Evidence = Gex.obs, DGex.obs',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.66*(left+right), 0.88*(bottom+top), 'meta.h posterior',horizontalalignment='left',verticalalignment='center',fontsize=18,fontname="Arial",fontweight="normal",color='red',transform=ax1.transAxes)
plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("meta_t2d_glp1.png", transparent=True,dpi=400)
