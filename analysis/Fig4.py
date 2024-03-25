# %%
import numpy as np
from os import listdir
from os.path import isfile, join
from matplotlib import ticker, cm
from matplotlib.pyplot import MultipleLocator
import matplotlib.pyplot as plt
from scipy.stats import norm
import glob
from analyze_metamodel import *


surrogate_ve_filepath = '../metamodel/Output/exocytosis_prior_wDGd.csv'
surrogate_pr_filepath = '../metamodel/Output/postprandial_prior_normal_wDGd_s1.csv'
surrogate_pa_filepath = '../metamodel/Output/pancreas_prior_wDGd.csv'
surrogate_ve_mean, surrogate_ve_std = get_surrogate_ve_model(surrogate_ve_filepath)
s_ve_en,s_pr_en,s_pa_en,surrogate_entropy = get_surrogate_entropy(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)

meta_filepath = '../metamodel/Output/fig4_test_Gpl/'
files_sorted = []
enumerate_mean = sorted(np.unique(np.array([item.split('/')[-1].split('_')[3] for item in glob.glob(meta_filepath+'metamodel_normal_w_*')])), key=float)
enumerate_cov = sorted(np.unique(np.array([item.split('/')[-1].split('_')[-1][:-4] for item in glob.glob(meta_filepath+'metamodel_normal_w_*')])), key=float)
for mean in enumerate_mean:
    for cov in enumerate_cov:    
        files_sorted += ['metamodel_normal_w_'+mean+'_cov_'+cov+'.csv']

Gpl_mean = np.unique(np.array([float(file.split('_')[3]) for file in files_sorted]))
Gpl_cov = np.unique(np.array([float(file.split('_')[-1].split('csv')[0][:-1]) for file in files_sorted]))
print(Gpl_mean)
print(Gpl_cov)

selected_time_span = 420
sve_idx = -1

meta_sve_mean, meta_sve_std = [],[]
for file in files_sorted:
    meta_ve_mean, meta_ve_std, _, _, _, _ = get_metamodel(meta_filepath+file)
    meta_sve_mean += [meta_ve_mean[:selected_time_span,sve_idx]]
    meta_sve_std += [meta_ve_std[:selected_time_span,sve_idx]]
meta_sve_mean = np.array(meta_sve_mean)
meta_sve_std = np.array(meta_sve_std)



# %% ##################################### Fig 4E #####################################

meta_entropy_list = np.recfromtxt('./fig4_test_Gpl_entropy.txt')

x = Gpl_mean*5.1 + (1-Gpl_mean[5])*5.1
y = np.linspace(0, 10, 11)
X, Y = np.meshgrid(x, y)
Z = []
for i in range(121):
    temp = np.array(meta_entropy_list)[i,:]
    # temp = np.array(meta_entropy_list)[i]
    temp = temp[~np.isnan(temp)]
    Z += [np.mean(np.array(temp))-np.mean(surrogate_entropy)]
Z = np.array(Z).reshape(11,11).transpose()


# %%
fig, ax = plt.subplots(figsize=(7, 6))

# levels = np.arange(np.min(Z), np.max(Z), (np.max(Z)-np.min(Z))/5000)
levels = np.arange(np.min(Z), np.max(Z), (np.max(Z)-np.min(Z))/5000)
cs = ax.contourf(X, Y, Z, levels, cmap=plt.get_cmap('coolwarm'))
cbar = fig.colorbar(cs, pad=0.01)
cbar.ax.ticklabel_format(style='sci',scilimits=(0,0),axis='x')
cbar.locator = ticker.MaxNLocator(nbins=6)
cbar.ax.tick_params(labelsize=25)
cbar.update_ticks()
# cbar.set_ticks(np.arange(round(np.min(Z),1), round(np.max(Z),1)+1, 1))
formatter = ticker.FormatStrFormatter('%.1f')
cbar.ax.yaxis.set_major_formatter(formatter)

ax.set_xlabel('<$p(G^{C}_{pl}|G^{PR})$>',fontsize = 25, fontname = 'Arial Narrow')
ax.set_ylabel('log$_{10}\ [\sigma^2$ of $p(G^{C}_{pl}|G^{PR})$]',fontsize = 25, 
              fontname = 'Arial Narrow', labelpad=0.1)

plt.yticks(np.arange(0, 11, 2), [round(i, 2) for i in np.linspace(-3, 4, 11)[::2]])
plt.xticks(np.arange(2.55, 8.16, 0.51*2)[:-1]+0.4, 
           [round(i, 1)+0.4 for i in np.linspace(2.55, 7.65, 11)[::2]][:-1])
ax.annotate(r'$\sum h(\mathbf{X^{i}_{s,upd}}) - \sum h(\mathbf{X^{i}_{s}})$', 
             xy=(0.1, 1.05), xycoords="axes fraction", fontsize=25)
ax.tick_params(labelsize=25)

plt.savefig('./Fig4_H.png', dpi=600, bbox_inches = 'tight')
# plt.close()





# %% ##################################### Fig 4F #####################################

# f = open('./fig4_test_Gpl_consistency.txt', 'w')
# for file in files_sorted:
#     model_consistency = get_model_consistency(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,meta_filepath+file)
#     print(file, np.mean(model_consistency), file=f)
# f.close()

meta_consistency = []
f = open('./fig4_test_Gpl_consistency.txt').readlines()
for line in f:
    if '.csv' in line:
        meta_consistency += [float(line.split(' ')[1].split('\n')[0])]



x = Gpl_mean*5.1 + (1-Gpl_mean[5])*5.1
y = np.linspace(0, 10, 11)
X, Y = np.meshgrid(x, y)
Z = np.array(meta_consistency).reshape(11,11).transpose()
    
fig, ax = plt.subplots(figsize=(7, 6))

# levels = np.arange(np.min(Z), np.max(Z), (np.max(Z)-np.min(Z))/5000)
levels = np.arange(0.36, 0.61, (np.max(Z)-np.min(Z))/5000)
cs = ax.contourf(X, Y, Z, levels, cmap=plt.get_cmap('coolwarm'))
cbar = fig.colorbar(cs, pad=0.01)
cbar.ax.ticklabel_format(style='sci',scilimits=(0,0),axis='x')
cbar.locator = ticker.MaxNLocator(nbins=6)
cbar.ax.tick_params(labelsize=25)
cbar.update_ticks()

ax.set_xlabel('<$p(G^{C}_{pl}|G^{PR})$>',fontsize = 25, fontname = 'Arial Narrow')
ax.set_ylabel('log$_{10}\ [\sigma^2$ of $p(G^{C}_{pl}|G^{PR})$]',fontsize = 25, fontname = 'Arial Narrow', labelpad=0.1)

plt.yticks(np.arange(0, 11, 2), [round(i, 2) for i in np.linspace(-3, 4, 11)[::2]])
plt.xticks(np.arange(2.55, 8.16, 0.51*2)[:-1]+0.4, 
           [round(i, 1)+0.4 for i in np.linspace(2.55, 7.65, 11)[::2]][:-1])
# ax.set_title('Model consistency between surrogate\n and updated surrogate model', fontsize=25)
ax.annotate(r'$\bar{\eta}$' + r'($\mathbf{X_{s,upd}}, \mathbf{X_{s}})$', 
            xy=(0.2, 1.05), xycoords="axes fraction", fontsize=25)
ax.tick_params(labelsize=25)

plt.savefig('./Fig4_I.png', dpi=200, bbox_inches = 'tight' )
plt.show()





# %% ##################################### Fig 4H #####################################

meta_file_H = 'metamodel_normal_w_0.2_cov_0.63096.csv'
meta_ve_mean, meta_ve_std, _, _, _, _ = get_metamodel(meta_filepath+meta_file_H)

timestep = 80
meta_sve_mean_H = meta_ve_mean[timestep, sve_idx]
meta_sve_std_H = meta_ve_std[timestep, sve_idx]
surrogate_sve_mean_H = surrogate_ve_mean[timestep, sve_idx]
surrogate_sve_std_H = surrogate_ve_std[timestep, sve_idx]

fig =plt.figure(figsize=(4, 3))

plot_normpdf(u=surrogate_sve_mean_H, sig=surrogate_sve_std_H, 
             color='C2', style='--',xlabel='$S^{VE}$',ylim=(-3e+7, 1e+9), 
             label='surrogate model')
plot_normpdf(u=meta_sve_mean_H, sig=meta_sve_std_H, 
             color='indianred', label='updated surrogate model')
plt.annotate('$t$=80min',xy=(0.05,0.9), xycoords="axes fraction", fontsize=25)
print(compute_overlap_of_normal_dist(surrogate_sve_mean_H, meta_sve_mean_H, surrogate_sve_std_H, meta_sve_std_H))
plt.annotate('0.75',xy=(0.4,0.1), xycoords="axes fraction", fontsize=23)

plt.subplots_adjust(wspace=0.13)
# plt.legend(prop={'size': 22, 'family': 'Arial Narrow'}, frameon=False)
plt.savefig('./Fig4_H.png',dpi=600, bbox_inches = 'tight')
plt.show()

# %%
timestep = 300
meta_sve_mean_H = meta_ve_mean[timestep, sve_idx]
meta_sve_std_H = meta_ve_std[timestep, sve_idx]
surrogate_sve_mean_H = surrogate_ve_mean[timestep, sve_idx]
surrogate_sve_std_H = surrogate_ve_std[timestep, sve_idx]

fig =plt.figure(figsize=(4, 3))

plot_normpdf(u=surrogate_sve_mean_H, sig=surrogate_sve_std_H, 
             color='C2', style='--',xlabel='$S^{VE}$',ylim=(-3e+7, 1e+9), 
             label='Surrogate model')
plot_normpdf(u=meta_sve_mean_H, sig=meta_sve_std_H, 
             color='indianred', label='Updated surrogate model')
plt.annotate('$t$=300min',xy=(0.05,0.9), xycoords="axes fraction", fontsize=25)
print(compute_overlap_of_normal_dist(surrogate_sve_mean_H, meta_sve_mean_H, surrogate_sve_std_H, meta_sve_std_H))
plt.annotate('0.77',xy=(0.4,0.1), xycoords="axes fraction", fontsize=23)

plt.subplots_adjust(wspace=0.13)
plt.legend(loc=(1,0), prop={'size': 22, 'family': 'Arial Narrow'}, frameon=False)
plt.savefig('./Fig4_H_300.png',dpi=600, bbox_inches = 'tight')
plt.show()
# %%
