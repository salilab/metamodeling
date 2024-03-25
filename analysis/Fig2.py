# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker 
import matplotlib as mpl
import numpy as np
import scipy.stats
from matplotlib.pyplot import MultipleLocator
import math
from analyze_metamodel import *
from scipy import stats

############################# obtain the joint distribution of input model ###########################################

I_pl_exp_data = './input_pr/plasma_insulin_obs.dat'
plasma_insulin = np.array([(float(item.split('\t')[1]), float(item.split('\t')[2])) for item in open(I_pl_exp_data).read().split('\n')[1:]])
plasma_insulin_std_ratio = np.mean(plasma_insulin[:, 1] / plasma_insulin[:, 0])

S_exp_data = './input_pr/insulin_secretion_obs.dat'
insulin_secretion = np.array([(float(item.split('\t')[1]), float(item.split('\t')[2]))  for item in open(S_exp_data).read().split('\n')[1:]])
insulin_secretion_std_ratio = np.mean(insulin_secretion[:, 1] / insulin_secretion[:, 0])

G_exp_data = './input_pr/plasma_glucose_obs.dat'
plasma_glucose = np.array([(float(item.split('\t')[1]), float(item.split('\t')[2]))  for item in open(G_exp_data).read().split('\n')[1:]])
plasma_glucose_std_ratio = np.mean(plasma_glucose[:, 1] / plasma_glucose[:, 0])

pr_std_avg_ratio = np.mean([plasma_insulin_std_ratio, insulin_secretion_std_ratio, plasma_glucose_std_ratio])
pr_input_mean = np.genfromtxt('./input_pr/input_pr_12_var_mean.txt')
pr_input_std = abs(pr_input_mean)*pr_std_avg_ratio
std_ratio = np.mean([plasma_insulin_std_ratio, insulin_secretion_std_ratio, plasma_glucose_std_ratio])


pr_input_std[:, 6] = abs(pr_input_mean[:, 6])*plasma_insulin_std_ratio
pr_input_std[:, 4] = abs(pr_input_mean[:, 4])*insulin_secretion_std_ratio
pr_input_std[:, 2] = abs(pr_input_mean[:, 2])*plasma_glucose_std_ratio

input_pr_entropy = compute_model_entropy(pr_input_mean, pr_input_std)
print('input model entropy: ', np.mean(input_pr_entropy[~np.isnan(input_pr_entropy)]))

input_8var_pr_entropy = compute_model_entropy(pr_input_mean[:, :8], pr_input_std[:, :8])
print('input model 8 var entropy: ', np.mean(input_8var_pr_entropy[~np.isnan(input_8var_pr_entropy)]))


# %% ############################# get model entropy and model consistency ###########################################

surrogate_pr_1 = '../metamodel/Output/postprandial_prior_normal_v1.csv'
surrogate_pr_mean_1, surrogate_pr_std_1 = get_surrogate_pr_model(surrogate_pr_1)
pr_entropy1 = compute_model_entropy(surrogate_pr_mean_1, surrogate_pr_std_1)
print('surrogate model 1 entropy: ', np.mean(pr_entropy1))

# %%

surrogate_pr_2 = '../metamodel/Output/postprandial_prior_normal_v2.csv'
surrogate_pr_mean_2, surrogate_pr_std_2 = get_surrogate_pr_model(surrogate_pr_2)
pr_entropy2 = compute_model_entropy(surrogate_pr_mean_2, surrogate_pr_std_2)
print('surrogate model 2 entropy: ', np.mean(pr_entropy2))

# %%

pr_consistency1, pr_consistency2 = [], []

plasma_insulin_idx = -2
insulin_secretion_idx = 4
plasma_glucose_idx = 2

input_G_idx = 2
input_S_idx = 4
input_I_idx = 6

for ts in range(len(surrogate_pr_mean_1)):

    surroagte_pr_I_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,plasma_insulin_idx], pr_input_mean[ts, 3], 
                                           surrogate_pr_std_1[ts,plasma_insulin_idx], pr_input_std[ts, 3])
    surroagte_pr_S_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,insulin_secretion_idx]/78, pr_input_mean[ts, 10], 
                                           surrogate_pr_std_1[ts,insulin_secretion_idx]/78, pr_input_mean[ts, 10])
    surroagte_pr_G_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, 0], 
                                           surrogate_pr_std_1[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, 0])
    pr_consistency1 += [np.mean([surroagte_pr_I_consistency1, surroagte_pr_S_consistency1, surroagte_pr_G_consistency1])]

    surroagte_pr_I_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,plasma_insulin_idx], pr_input_mean[ts, 3], 
                                           surrogate_pr_std_2[ts,plasma_insulin_idx], pr_input_std[ts, 3])
    surroagte_pr_S_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,insulin_secretion_idx]/78, pr_input_mean[ts, 10], 
                                           surrogate_pr_std_2[ts,insulin_secretion_idx]/78, pr_input_mean[ts, 10])
    surroagte_pr_G_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, 0], 
                                           surrogate_pr_std_2[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, 0])
    pr_consistency2 += [np.mean([surroagte_pr_I_consistency2, surroagte_pr_S_consistency2, surroagte_pr_G_consistency2])]

print('surrogate model 1 consistency: ', np.mean(pr_consistency1))
print('surrogate model 2 consistency: ', np.mean(pr_consistency2))



pr_consistency1, pr_consistency2 = [], []

plasma_insulin_idx = -2
insulin_secretion_idx = 4
plasma_glucose_idx = 2

for ts in range(len(surrogate_pr_mean_1)):

    surroagte_pr_I_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,plasma_insulin_idx], pr_input_mean[ts, input_I_idx], 
                                           surrogate_pr_std_1[ts,plasma_insulin_idx], pr_input_std[ts, input_I_idx])
    surroagte_pr_S_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,insulin_secretion_idx]/78, pr_input_mean[ts, input_S_idx], 
                                           surrogate_pr_std_1[ts,insulin_secretion_idx]/78, pr_input_std[ts, input_S_idx])
    surroagte_pr_G_consistency1 = compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, input_G_idx], 
                                           surrogate_pr_std_1[ts,plasma_glucose_idx]*18.016, pr_input_std[ts, input_G_idx])
    
    surroagte_pr_var_consistency1 = [compute_overlap_of_normal_dist(surrogate_pr_mean_1[ts, i], pr_input_mean[ts, i], 
                                                                    surrogate_pr_std_1[ts, i], pr_input_std[ts, i]) for i in range(8)]

    surroagte_pr_var_consistency1[insulin_secretion_idx] = surroagte_pr_S_consistency1
    surroagte_pr_var_consistency1[plasma_glucose_idx] = surroagte_pr_G_consistency1

    pr_consistency1 += [np.mean(surroagte_pr_var_consistency1)]

    surroagte_pr_I_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,plasma_insulin_idx], pr_input_mean[ts, input_I_idx], 
                                           surrogate_pr_std_2[ts,plasma_insulin_idx], pr_input_std[ts, input_I_idx])
    surroagte_pr_S_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,insulin_secretion_idx]/78, pr_input_mean[ts, input_S_idx], 
                                           surrogate_pr_std_2[ts,insulin_secretion_idx]/78, pr_input_std[ts, input_S_idx])
    surroagte_pr_G_consistency2 = compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts,plasma_glucose_idx]*18.016, pr_input_mean[ts, input_G_idx], 
                                           surrogate_pr_std_2[ts,plasma_glucose_idx]*18.016, pr_input_std[ts, input_G_idx])
    
    surroagte_pr_var_consistency2 = [compute_overlap_of_normal_dist(surrogate_pr_mean_2[ts, i], pr_input_mean[ts, i], 
                                                                    surrogate_pr_std_2[ts, i], pr_input_std[ts, i]) for i in range(8)]

    surroagte_pr_var_consistency2[insulin_secretion_idx] = surroagte_pr_S_consistency2
    surroagte_pr_var_consistency2[plasma_glucose_idx] = surroagte_pr_G_consistency2

    pr_consistency2 += [np.mean(surroagte_pr_var_consistency2)]

print('surrogate model 1 consistency: ', np.mean(pr_consistency1))
print('surrogate model 2 consistency: ', np.mean(pr_consistency2))








# %% ########################## Fig 2B ##############################################

fig, ax= plt.subplots(1,1, figsize=(9,4))

time_data = [0,   5,  10,  20,  30,  40,  50,  60,  75,  90, 105, 120, 150, 180, 210, 240, 260, 280, 300, 360, 420]

xrange = np.arange(0,420,1)

x_major_locator=MultipleLocator(100)
y_major_locator=MultipleLocator(200)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
x_minor_locator=MultipleLocator(50)
y_minor_locator=MultipleLocator(100)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.set_ylim(-101,650)
ax.set_xlabel('Time [min]', fontproperties = 'Arial Narrow', size = 23)
ax.set_ylabel('$I^{PR}$ [pM]', fontproperties = 'Arial Narrow', size = 23)
# ax.plot(xrange, pr_input_mean[:, 3], 'grey', label='Input model, \n $N$=12, $h$=115.2', 
#            linewidth=2, linestyle = "dashed")
input_pr_I_entropy = compute_model_entropy(pr_input_mean[:, input_I_idx].reshape(-1,1), 
                                           pr_input_std[:, input_I_idx].reshape(-1,1))
ax.plot(xrange, pr_input_mean[:, input_I_idx], 'grey', label=r'$h(I^{PR}_{\text{in}})$='+'{:.1f}'.format(np.mean(input_pr_I_entropy)), 
           linewidth=2, linestyle = "dashed")
ax.fill_between(xrange, pr_input_mean[:, input_I_idx]-pr_input_std[:, input_I_idx], 
                    pr_input_mean[:, input_I_idx]+pr_input_std[:, input_I_idx], alpha=0.2, color='grey')
surroagte_pr_I_entropy1 = compute_model_entropy(surrogate_pr_mean_1[:, plasma_insulin_idx].reshape(-1,1), 
                                                surrogate_pr_std_1[:, plasma_insulin_idx].reshape(-1,1))

ax.plot(xrange, surrogate_pr_mean_1[:, plasma_insulin_idx], 
           label=r'$h(I_{1,s}^{PR})$='+'{:.1f}, '.format(np.mean(surroagte_pr_I_entropy1)) + \
            r'$\eta(I_{1,s}^{PR}, I_{\text{in}}^{PR})$='+'{:.2f}'.format(np.mean(surroagte_pr_I_consistency1)), 
           color='C2', linewidth=2)
ax.fill_between(xrange, surrogate_pr_mean_1[:, plasma_insulin_idx]-surrogate_pr_std_1[:, plasma_insulin_idx], 
                          surrogate_pr_mean_1[:, plasma_insulin_idx]+surrogate_pr_std_1[:, plasma_insulin_idx],
                  color='C2',  alpha=0.2)
surroagte_pr_I_entropy2 = compute_model_entropy(surrogate_pr_mean_2[:, plasma_insulin_idx].reshape(-1,1), 
                                                surrogate_pr_std_2[:, plasma_insulin_idx].reshape(-1,1))
ax.plot(xrange, surrogate_pr_mean_2[:, plasma_insulin_idx], 
           label=r'$h(I_{2,s}^{PR})$='+'{:.1f}, '.format(np.mean(surroagte_pr_I_entropy2)) + \
            r'$\eta(I_{2,s}^{PR}, I_{\text{in}}^{PR})$='+'{:.2f}'.format(np.mean(surroagte_pr_I_consistency2)), 
           color='C0', linewidth=2)
ax.fill_between(xrange, surrogate_pr_mean_2[:, plasma_insulin_idx]-surrogate_pr_std_2[:, plasma_insulin_idx], 
                          surrogate_pr_mean_2[:, plasma_insulin_idx]+surrogate_pr_std_2[:, plasma_insulin_idx],
                color='C0',  alpha=0.2)
ax.legend(loc='upper right', prop={'size': 20, 'family': 'Arial Narrow'}, frameon=False)
ax.vlines(30, -100, surrogate_pr_mean_2[30, plasma_insulin_idx],linestyles = "dashed",color='darkred')

plt.savefig('./Fig2_B.png',dpi=600, bbox_inches = 'tight')
plt.show()

# %%

fig, ax= plt.subplots(1,1, figsize=(9,4))

ax.set_xlabel('Time [min]', fontproperties = 'Arial Narrow', size = 23)
ax.set_ylabel('$S^{PR}$ [pmol/kg/min]', fontproperties = 'Arial Narrow', size = 23, labelpad=10)
x_major_locator=MultipleLocator(100)
y_major_locator=MultipleLocator(5)
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
x_minor_locator=MultipleLocator(50)
y_minor_locator=MultipleLocator(2.5)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.set_ylim(-1,19)
input_pr_S_entropy = compute_model_entropy(pr_input_mean[:, input_S_idx].reshape(-1,1), 
                                           pr_input_std[:, input_S_idx].reshape(-1,1))
ax.plot(xrange, pr_input_mean[:, input_S_idx], 'grey', label=r'$h(S_{in}^{PR})$'+'={:.1f}'.format(np.mean(input_pr_S_entropy)), 
           linewidth=2, linestyle = "dashed")
ax.fill_between(xrange, pr_input_mean[:, input_S_idx]-pr_input_std[:, input_S_idx], 
                    pr_input_mean[:, input_S_idx]+pr_input_std[:, input_S_idx], alpha=0.2, color='grey')

surroagte_pr_S_entropy1 = compute_model_entropy(surrogate_pr_mean_1[:, insulin_secretion_idx].reshape(-1,1), 
                                                surrogate_pr_std_1[:, insulin_secretion_idx].reshape(-1,1))
ax.plot(xrange, surrogate_pr_mean_1[:, insulin_secretion_idx]/78, 'C2', 
           label=r'$h(S_{1,s}^{PR})$='+'{:.1f}, '.format(np.mean(surroagte_pr_S_entropy1)) + \
            r'$\eta(S_{1,s}^{PR}, S_{\text{in}}^{PR})$='+'{:.2f}'.format(np.mean(surroagte_pr_S_consistency1)), 
           linewidth=2)
ax.fill_between(xrange, surrogate_pr_mean_1[:, insulin_secretion_idx]/78-surrogate_pr_std_1[:, insulin_secretion_idx]/78, 
                        surrogate_pr_mean_1[:, insulin_secretion_idx]/78+surrogate_pr_std_1[:, insulin_secretion_idx]/78,
                color='C2',  alpha=0.2)

surroagte_pr_S_entropy2 = compute_model_entropy(surrogate_pr_mean_2[:, insulin_secretion_idx].reshape(-1,1), 
                                                surrogate_pr_std_2[:, insulin_secretion_idx].reshape(-1,1))
ax.plot(xrange, surrogate_pr_mean_2[:, insulin_secretion_idx]/78, 'C0', 
           label=r'$h(S_{2,s}^{PR})$='+'{:.1f}, '.format(np.mean(surroagte_pr_S_entropy2)) + \
            r'$\eta(S_{2,s}^{PR}, S_{\text{in}}^{PR})$='+'{:.2f}'.format(np.mean(surroagte_pr_S_consistency2)), 
           linewidth=2)
ax.fill_between(xrange, surrogate_pr_mean_2[:, insulin_secretion_idx]/78-surrogate_pr_std_2[:, insulin_secretion_idx]/78, 
                        surrogate_pr_mean_2[:, insulin_secretion_idx]/78+surrogate_pr_std_2[:, insulin_secretion_idx]/78,
                color='C0',   alpha=0.2)
ax.legend(loc='upper right', prop={'size': 20, 'family': 'Arial Narrow'}, frameon=False)
ax.vlines(300, -1, surrogate_pr_mean_1[300, insulin_secretion_idx]/78, linestyles = "dashed",color='darkred')


plt.savefig('./Fig2_E.png',dpi=600, bbox_inches = 'tight')
plt.show()




# %% ########################## Fig 2CD ##############################################

fig =plt.figure(figsize=(4, 9))

ax = fig.add_subplot(2, 1, 1)

timestep=30

plot_normpdf(u=pr_input_mean[timestep, input_I_idx], sig=plasma_insulin[4, 1], 
             color='grey', xlabel='$I^{PR}$', ylim=(-0.001, 0.015), style='--')
plot_normpdf(u=surrogate_pr_mean_1[timestep, plasma_insulin_idx], 
             sig=surrogate_pr_std_1[timestep, plasma_insulin_idx], 
             color='C2')
v1 = compute_overlap_of_normal_dist(pr_input_mean[timestep, input_I_idx], surrogate_pr_mean_1[timestep, plasma_insulin_idx], 
                                     pr_input_std[timestep, input_I_idx], surrogate_pr_std_1[timestep, plasma_insulin_idx])
print(v1)
plot_normpdf(u=surrogate_pr_mean_2[timestep, plasma_insulin_idx], 
             sig=surrogate_pr_std_2[timestep, plasma_insulin_idx], 
             color='C0')
v2 = compute_overlap_of_normal_dist(pr_input_mean[timestep, input_I_idx], surrogate_pr_mean_2[timestep, plasma_insulin_idx], 
                                     pr_input_std[timestep, input_I_idx], surrogate_pr_std_2[timestep, plasma_insulin_idx])
print(v2)
ax.annotate('$\eta_1$={:.2f}'.format(v1),xy=(0.2,0.7), xycoords="axes fraction", fontsize=22,color='C2')
ax.annotate('$\eta_2$={:.2f}'.format(v2),xy=(0.65,0.55), xycoords="axes fraction", fontsize=22,color='C0')
ax.annotate('$t$=30min, $\sigma_1$={:.1f}, $\sigma_2$={:.1f}'.format(surrogate_pr_std_1[timestep, 
                                                                                        plasma_insulin_idx], 
                                                                     surrogate_pr_std_2[timestep, 
                                                                                        plasma_insulin_idx]),
                                                                     xy=(0.05,0.9), xycoords="axes fraction", 
                                                                     fontsize=22)



ax = fig.add_subplot(2, 1, 2)

timestep=300

plot_normpdf(u=pr_input_mean[timestep, input_S_idx], sig=pr_input_std[timestep, input_S_idx], 
             color='grey', style='--', ylim=(-0.05, 0.72), xlabel='$S^{PR}$')
plot_normpdf(u=surrogate_pr_mean_1[timestep, insulin_secretion_idx]/78, 
             sig=surrogate_pr_std_1[timestep, insulin_secretion_idx]/78, 
             color='C2')
v1 = compute_overlap_of_normal_dist(pr_input_mean[timestep, input_S_idx], surrogate_pr_mean_1[timestep, insulin_secretion_idx]/78, 
                                     pr_input_std[timestep, input_S_idx], surrogate_pr_std_1[timestep, insulin_secretion_idx]/78)
print(v1)
plot_normpdf(u=surrogate_pr_mean_2[timestep, insulin_secretion_idx]/78, 
             sig=surrogate_pr_std_2[timestep, insulin_secretion_idx]/78, 
             color='C0')
v2 = compute_overlap_of_normal_dist(pr_input_mean[timestep, input_S_idx], surrogate_pr_mean_2[timestep, insulin_secretion_idx]/78, 
                                     pr_input_std[timestep, input_S_idx], surrogate_pr_std_2[timestep, insulin_secretion_idx]/78)
print(v2)
ax.annotate('$\eta_1$={:.2f}'.format(v1),xy=(0.65,0.6), xycoords="axes fraction", fontsize=22,color='C2')
ax.annotate('$\eta_2$={:.2f}'.format(v2),xy=(0.45,0.75), xycoords="axes fraction", fontsize=22,color='C0')
ax.annotate('$t$=300min, $\sigma_1$={:.1f}, $\sigma_2$={:.1f}'.format(surrogate_pr_std_1[timestep, insulin_secretion_idx]/78,
                                                                      surrogate_pr_std_2[timestep, insulin_secretion_idx]/78),
                                                                      xy=(0.05,0.9), xycoords="axes fraction", fontsize=22)

plt.subplots_adjust(wspace=0.1)
plt.savefig('./Fig2_CD.png',dpi=600, bbox_inches = 'tight' )

# %%
