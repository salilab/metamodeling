# %%
import numpy as np
from scipy.stats import *
import matplotlib.pyplot as plt
from scipy.stats import entropy
import statistics_basic as stat
import matplotlib.pyplot as plt
import matplotlib.ticker 
import matplotlib as mpl
import numpy as np
from matplotlib.pyplot import MultipleLocator
from analyze_metamodel import *
import glob

entropy_hist = []
f = open('./fig3_full_time_entropy.txt').readlines()
for line in f:
    if '.csv' in line:
        entropy_hist += [float(line.split(' ')[1].split('\n')[0])]


surrogate_ve_filepath = '../metamodel/Output/exocytosis_prior_wDGd.csv'
surrogate_pr_filepath = '../metamodel/Output/postprandial_prior_normal_wDGd.csv'
surrogate_pa_filepath = '../metamodel/Output/pancreas_prior_wDGd.csv'
s_ve_en,s_pr_en,s_pa_en,surrogate_entropy = get_surrogate_entropy(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)


# %%  ##################################### Fig 3B #####################################

fig, ax= plt.subplots(figsize=(5,5))

mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20) 
ax.set_xlabel('Model entropy', fontproperties = 'Arial Narrow', size = 25)
ax.set_ylim(0,0.26)
ax.tick_params(labelsize=25)
# ax.set_xlim(175.1,179.1)
x_major_locator=MultipleLocator(1)
y_major_locator=MultipleLocator(0.05)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.f'))
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
# x_minor_locator=MultipleLocator(0.1)
y_minor_locator=MultipleLocator(0.025)
# ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_minor_locator(y_minor_locator)

mean_val = np.nanmean(entropy_hist)  # Calculate the mean, ignoring NaN
entropy_hist = np.where(np.isnan(entropy_hist), mean_val, entropy_hist) 

n, bins = np.histogram(np.array(entropy_hist), bins=20)
n_edited = n / len(entropy_hist)  
bin_width = np.diff(bins)[0]
bin_centers = (bins[:-1] + bins[1:]) / 2
bars = plt.bar(bin_centers, n_edited, width=bin_width, color='gray', label='After metamodeling', alpha=0.6)
plt.axvline(np.mean(surrogate_entropy), 0, 0.66, linestyle="dashed", linewidth=2, color='darkred', label='Before metamodeling')
plt.annotate(str(round(np.mean(surrogate_entropy),1)), color='darkred', xy=(0.06, 0.68), xycoords="axes fraction", fontsize=22)

for bar in bars:
    plt.gca().add_patch(plt.Rectangle((bar.get_x(), bar.get_y()), bar.get_width(), bar.get_height(), fill=False, edgecolor='black', lw=1))


plt.legend(loc='upper center', prop={'size': 22, 'family': 'Arial Narrow'}, frameon=False)
plt.savefig('./Fig3_B.png', dpi=600, bbox_inches='tight')
plt.show()

# %%

meta_files = glob.glob('../metamodel/Output/fig3_full_time_CPDs/*')
f = open('./fig3_full_time_consistency.txt', 'w')

for file in meta_files:
    model_consistency = get_model_consistency(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,file)
    print(file, np.mean(model_consistency), file=f)

f.close()

# %%  ##################################### Fig 3C #####################################

consistency_hist = []
f = open('./fig3_full_time_consistency.txt').readlines()
for line in f:
    if '.csv' in line and '_0_' not in line:
        consistency_hist += [float(line.split(' ')[1].split('\n')[0])]

fig, ax= plt.subplots(figsize=(5,5))

mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20) 
ax.set_xlabel('Model consistency', fontproperties = 'Arial Narrow', size = 25)
# ax.set_title('Model consistency distribution for all \n surrogate models using different CPDs', fontproperties = 'Arial Narrow', size = 23)
# ax.set_ylabel('Frequency', fontproperties = 'Arial Narrow', size = 25)
# ax.set_title('Metamodel', fontproperties = 'Arial Narrow', size = 30)
ax.set_ylim(0,0.21)
ax.tick_params(labelsize=25)
ax.set_xlim(0.44,0.66)
x_major_locator=MultipleLocator(0.05)
y_major_locator=MultipleLocator(0.05)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
x_minor_locator=MultipleLocator(0.025)
ax.xaxis.set_minor_locator(x_minor_locator)

n, bins = np.histogram(np.array(consistency_hist), bins=20)
n_edited = n / len(consistency_hist)  
bin_width = np.diff(bins)[0]
bin_centers = (bins[:-1] + bins[1:]) / 2
bars = plt.bar(bin_centers, n_edited, width=bin_width, color='grey', label='Metamodel', alpha=0.6)
for bar in bars:
    plt.gca().add_patch(plt.Rectangle((bar.get_x(), bar.get_y()), bar.get_width(), bar.get_height(), fill=False, edgecolor='black', lw=1))

# plt.legend(loc='upper center', prop={'size': 20, 'family': 'Arial Narrow'}, frameon=False)
plt.savefig('./Fig3_C.png', dpi=600, bbox_inches='tight')
plt.show()



# %%  ##################################### Fig 3D #####################################


consistency_hist = []
f = open('./fig3_full_time_consistency.txt').readlines()
for line in f:
    if '.csv' in line:
        consistency_hist += [(line.split(' ')[0].split('/')[-1] , float(line.split(' ')[1].split('\n')[0]))]

entropy_hist = []
f = open('./fig3_full_time_entropy.txt').readlines()
for line in f:
    if '.csv' in line:
        entropy_hist += [(line.split(' ')[0].split('/')[-1] , float(line.split(' ')[1].split('\n')[0]))]

consistency_, entropy_ = [], []
consistency_hist_ = [item[0] for item in consistency_hist]
for item in entropy_hist:
    filename = item[0]
    if '_0_' not in filename:
        entropy_ += [item[1]]
        consistency_ += [consistency_hist[consistency_hist_.index(filename)][1]]


fig, ax= plt.subplots(figsize=(5,5))

mpl.rc('xtick', labelsize=25) 
mpl.rc('ytick', labelsize=25) 
ax.tick_params(labelsize=25)
ax.set_xlabel('Model entropy', fontproperties = 'Arial Narrow', size = 23)
# ax.set_title('Change of model consistency as \n the change of model entropy', fontproperties = 'Arial Narrow', size = 23)
ax.set_ylabel('Model consistency', fontproperties = 'Arial Narrow', size = 23)
# ax.set_title('Metamodel', fontproperties = 'Arial Narrow', size = 30)
ax.set_ylim(0.44,0.66)
# ax.set_xlim(0.45,0.49)
x_major_locator=MultipleLocator(1)
y_major_locator=MultipleLocator(0.05)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.f'))
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
# x_minor_locator=MultipleLocator(0.1)
# ax.xaxis.set_minor_locator(x_minor_locator)

plt.scatter(entropy_, consistency_, s=14, color='gray', alpha=0.6)
plt.axvline(np.mean(surrogate_entropy), 0, 1, linestyle="dashed", linewidth=2, color='darkred', label='Before metamodeling')
# plt.legend(loc='upper center', prop={'size': 20, 'family': 'Arial Narrow'}, frameon=False)
plt.savefig('./Fig3_D.png', dpi=600, bbox_inches='tight')
plt.show()


# %%
