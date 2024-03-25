# %%
from matplotlib import ticker, cm
from matplotlib.pyplot import MultipleLocator
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib as mpl
import math
from os import listdir
from os.path import isfile, join
import statistics_basic as stat
from analyze_metamodel import *
import glob


filepath = '../metamodel/Output/fig5_test_Gpl/'
enumerate_mean = sorted(np.unique(np.array([item.split('/')[-1].split('_')[8] for item in glob.glob(filepath+'metamodel_scan_input_G_PR_test_mean_scale_*_cov_*')])), key=float)
enumerate_cov = sorted(np.unique(np.array([item.split('/')[-1].split('_')[-1][:-4] for item in glob.glob(filepath+'metamodel_scan_input_G_PR_test_mean_scale_*_cov_*')])), key=float)
print(enumerate_mean)
print(enumerate_cov)
files_sorted = []
for mean in enumerate_mean:
    for cov in enumerate_cov:    
        files_sorted += ['metamodel_scan_input_G_PR_test_mean_scale_'+mean+'_cov_'+cov+'.csv']

Gpl_mean = np.unique(np.array([float(item) for item in enumerate_mean]))
Gpl_cov = np.unique(np.array([float(item) for item in enumerate_cov]))

entropy_list = []
f = open('./fig5_test_Gpl_full_time_entropy.txt').readlines()
for file in files_sorted:
    for line in f:
        if file in line:
            print(file)
            entropy_list += [float(line.split(' ')[1].split('\n')[0])]
        
x = np.linspace(-2.55, 2.55, 11)
y = np.linspace(0, 10, 11)
X, Y = np.meshgrid(x, y)

selected_timespan = 420

surrogate_ve_filepath = '../metamodel/Output/exocytosis_prior_wDGd.csv'
surrogate_pr_filepath = '../metamodel/Output/postprandial_prior_normal_wDGd_s1.csv'
surrogate_pa_filepath = '../metamodel/Output/pancreas_prior_wDGd.csv'
s_ve_en,s_pr_en,s_pa_en,surrogate_entropy = get_surrogate_entropy(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)

# %%
Z = np.array(entropy_list).reshape(11,11).transpose()-np.mean(surrogate_entropy)
Z[np.isnan(Z)] = 0

# %% ##################################### Fig 5B #####################################


fig, ax = plt.subplots(figsize=(6, 5))
levels = np.arange(np.min(Z), np.max(Z), (np.max(Z)-np.min(Z))/2000)
cs = ax.contourf(X, Y, Z, levels, cmap=plt.get_cmap('coolwarm')) # cmap=plt.get_cmap('coolwarm_r'))
cbar = fig.colorbar(cs, pad=0.02)
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.ax.tick_params(labelsize=20)
cbar.update_ticks()

plt.yticks(np.linspace(0,10,6), [round(i,2) for n,i in enumerate(np.linspace(-4, 2, 11)) if n % 2 ==0])
plt.xticks(np.arange(-2.55, 2.55+0.51, 0.51*2)[:-1]+0.5, 
           [round(i+0.5, 1) for i in np.arange(-2.55, 2.55+0.51, 0.51)[::2]][:-1])
ax = plt.gca()
ax.tick_params(labelsize=20)
ax.set_xlabel('$\widebar{err}(G^{PR})$ [mM]',fontsize = 18, fontname = 'Arial Narrow')
ax.set_ylabel('log$_{10}$ [$\sigma^2 (G^{PR})$] [mM]',fontsize = 18, fontname = 'Arial Narrow')
ax.annotate(r'$\sum h(\mathbf{X^{i}_{s,upd}}) - \sum h(\mathbf{X^{i}_{s}})$',
            xy=(0.1, 1.05), xycoords="axes fraction", fontsize=22)

plt.savefig('./Fig5_B.png', dpi=600, bbox_inches = 'tight')
plt.show()



# %% ##################################### Fig 5C #####################################

f = open('./fig5_test_Gpl_consistency.txt', 'w')
for file in files_sorted:
    model_consistency = get_model_consistency(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,filepath+file)
    print(file, np.mean(model_consistency), file=f)
f.close()

meta_consistency = []
f = open('./fig5_test_Gpl_consistency.txt').readlines()
for line in f:
    if '.csv' in line:
        meta_consistency += [float(line.split(' ')[1].split('\n')[0])]

print(np.array(meta_consistency).shape)
Z = np.array(meta_consistency).reshape(11,11).transpose()
print(np.min(Z), np.max(Z))

fig, ax = plt.subplots(figsize=(6, 5))

levels = np.arange(np.min(Z), np.max(Z)+0.01, (np.max(Z)-np.min(Z))/2000)
# levels = np.arange(0.45, 0.66, (np.max(Z)-np.min(Z))/2000)
cs = ax.contourf(X, Y, Z, levels, cmap=plt.get_cmap('coolwarm'))
cbar = fig.colorbar(cs, pad=0.02)
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.ax.tick_params(labelsize=20)
cbar.update_ticks()
# cbar.set_ticks(np.arange(round(np.min(Z),2), round(np.max(Z),2), 0.05))
formatter = ticker.FormatStrFormatter('%.2f')
cbar.ax.yaxis.set_major_formatter(formatter)

plt.yticks(np.linspace(0,10,6), [round(i,2) for n,i in enumerate(np.linspace(-4, 4, 11)) if n % 2 ==0])
plt.xticks(np.arange(-2.55, 2.55+0.51, 0.51*2)[:-1]+0.5, 
           [round(i+0.5, 1) for i in np.arange(-2.55, 2.55+0.51, 0.51)[::2]][:-1])
ax = plt.gca()
ax.tick_params(labelsize=20)
ax.set_xlabel('$\widebar{err}(G^{PR})$ [mM]',fontsize = 18, fontname = 'Arial Narrow')
ax.set_ylabel('log$_{10}$ [$\sigma^2 (G^{PR})$] [mM]',fontsize = 18, fontname = 'Arial Narrow')
ax.annotate(r'$\bar{\eta}$' + r'($\mathbf{X_{s,upd}}, \mathbf{X_{s}})$',
            xy=(0.25, 1.05), xycoords="axes fraction", fontsize=22)

plt.savefig('./Fig5_C.png', dpi=600, bbox_inches = 'tight')
plt.show()

# %%
