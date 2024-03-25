# %%
import numpy as np
from scipy.stats import *
import matplotlib.pyplot as plt
import math
from scipy.stats import entropy
import urllib.request
import os
import statistics_basic as stat
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker 
import matplotlib as mpl
import numpy as np
from matplotlib.pyplot import MultipleLocator
from analyze_metamodel import *
import glob

t2d_meta_file = '../metamodel/Output/metamodel_t2d_opt.csv'
t2d_meta_ve_mean, t2d_meta_ve_std, t2d_meta_pr_mean, t2d_meta_pr_std, t2d_meta_pa_mean, t2d_meta_pa_std = get_metamodel(t2d_meta_file)

non_opt_t2d_metafile = '../metamodel/Output/metamodel_t2d_non_opt.csv'
non_opt_t2d_meta_ve_mean, non_opt_t2d_meta_ve_std, non_opt_t2d_meta_pr_mean, non_opt_t2d_meta_pr_std, non_opt_t2d_meta_pa_mean, non_opt_t2d_meta_pa_std = get_metamodel(non_opt_t2d_metafile)

meta_file = '../metamodel/Output/metamodel_normal_opt.csv'
meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_metamodel(meta_file)

non_opt_metafile = '../metamodel/Output/metamodel_normal_non_opt.csv'
non_opt_meta_ve_mean, non_opt_meta_ve_std, non_opt_meta_pr_mean, non_opt_meta_pr_std, non_opt_meta_pa_mean, non_opt_meta_pa_std = get_metamodel(non_opt_metafile)


# %%

fig, ax= plt.subplots(1,1, figsize=(6,4))

mpl.rc('xtick', labelsize=23) 
mpl.rc('ytick', labelsize=23) 

time = np.arange(0, 420, 1)

# Gbasal I kt ve
idx = -2
ax.set_xlim([10,140])
ax.set_ylim(0,3)
ax.set_xlabel('Time [min]', fontproperties = 'Arial Narrow', size = 20)
ax.set_ylabel('Fold change of $I^{PR}_{normal}$ / $I^{PR}_{t2d}$', fontproperties = 'Arial Narrow', size = 22, labelpad=10)
x_major_locator=MultipleLocator(30)
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
# y_major_locator=MultipleLocator(250)
# ax.yaxis.set_major_locator(y_major_locator)
# ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
# y_minor_locator=MultipleLocator(100)

x_minor_locator=MultipleLocator(10)
ax.xaxis.set_minor_locator(x_minor_locator)

exp_time = [30, 60, 90, 120]
exp_t2d = [21.84684684684685, 24.024024024024026, 24.54954954954955, 25.375375375375377]
exp_normal = [44.66966966966967, 45.945945945945944, 43.768768768768766, 30.480480480480484]

exp_fc = np.array(exp_normal) / np.array(exp_t2d)
meta_pr_fc = [meta_pr_mean[i, idx] / t2d_meta_pr_mean[i, idx] for i in exp_time]
non_opt_meta_pr_fc = [non_opt_meta_pr_mean[i+1, idx] / non_opt_t2d_meta_pr_mean[i+1, idx] for i in exp_time]

ax.bar(np.array(exp_time)-5, exp_fc, width=5, color='Gray', label='Experimental measurements', linewidth=2, alpha=0.8)
ax.bar(np.array(exp_time), meta_pr_fc, width=5, color='indianred', label='Optimized metamodel', linewidth=2, alpha=0.8)
ax.bar(np.array(exp_time)+5, non_opt_meta_pr_fc, width=5, color='C0', label='Non-optimized metamodel', linewidth=2, alpha=0.8)

x = [30, 60, 90, 120]
plt.xticks([30, 60, 90, 120], x, rotation=0, fontsize = 20)

ax.legend(loc='upper right', prop={'size': 15, 'family': 'Arial Narrow'}, frameon=False)
# ax[0].set_title('Time courses of model variables', size = 25)

plt.subplots_adjust(hspace=0)
plt.savefig('./Fig6.png', dpi=600, bbox_inches = 'tight' )
plt.show()

# %%
