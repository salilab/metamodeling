# ##########################################################################
# Author: Chenxi Wang (wangchx@shanghaitech.edu.cn)
# Date: 2023-04
# ##########################################################################

import numpy as np
from scipy.stats import norm
import pandas as pd
import scipy.stats
import math
import numpy as np
from os import listdir
from os.path import isfile, join
import os
import urllib.request
import statistics_basic as stat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker 
import matplotlib as mpl
import numpy as np
from matplotlib.pyplot import MultipleLocator
from scipy import stats


font = {'family' : 'Arial narrow', 'size'   : 25}
COLOR = '#202020'
mpl.rc('font', **font)
mpl.rc('xtick', labelsize=25) 
mpl.rc('ytick', labelsize=25) 
mpl.rcParams['text.color'] = COLOR
mpl.rcParams['axes.labelcolor'] = COLOR
mpl.rcParams['xtick.color'] = COLOR
mpl.rcParams['ytick.color'] = COLOR


def entropy_1(X):
    probs = [np.mean(X == c) for c in set(X)]
    return np.sum(-p * np.log2(p) for p in probs)



# def tsallis_entropy(X, q):
#     """ Compute Tsallis entropy with parameter q.
#         X is a dataset and q is a scalar.
#     """
#     if q == 1:
#         return entropy_1(X)  # fallback to Shannon entropy
#     else:
#         probs = [np.mean(X == c) for c in set(X)]
#         return (1 - np.sum([p**q for p in probs])) / (q - 1)



def compute_model_entropy(mean, std):

    time = mean.shape[0]
    num_var = mean.shape[1]
    model_entropy = []

    for t in range(time):
        model_entropy_t = 0
        for var_idx in range(num_var):
            u = mean[t, var_idx]
            sig = std[t, var_idx]
            # model_entropy_t += norm(u, sig).entropy()
            
            x = np.linspace(u - 3*sig, u + 3*sig, 500)
            y = np.exp(-(x - u) ** 2 / (2 * sig ** 2)) / (math.sqrt(2*math.pi)*sig)
            model_entropy_t += entropy_1(y)
        
        model_entropy.append(model_entropy_t)
        
    return np.array(model_entropy)


def compute_overlap_of_normal_dist(m1,m2,std1,std2):

    if std1==0: std1 += 1e-10
    if std2==0: std2 += 1e-10
    
    N1 = stat.NormalDist(m1, std1)
    N2 = stat.NormalDist(m2, std2)

    return N1.overlap(N2)


def get_KS_p_value(mean1, mean2, std1, std2):
     
     x = np.random.normal(mean1, std1, 2000)
     y = np.random.normal(mean2, std2, 2000)
     
     ks_statistic, p_value = stats.ks_2samp(x, y)
     
#      print("K-S statistic:", ks_statistic)
#      print("p-value:", p_value)

     return p_value


def plot_normpdf(u, sig, color, xlabel=None, ylim=None, style='-', label=None):
    
    x = np.linspace(u - 3*sig, u + 3*sig, 100000)
    y = np.exp(-(x - u) ** 2 / (2 * sig ** 2)) / (math.sqrt(2*math.pi)*sig)

    if xlabel is not None:
        plt.xlabel(xlabel, labelpad=14)
    plt.ylabel('Distribution', labelpad=6)
    ax=plt.gca()  
    ax.set_xticks([])
    ax.set_yticks([])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    if label is not None:
        plt.plot(x, y, linestyle=style, linewidth=3, c=color, label=label)   
    else:
        plt.plot(x, y, linestyle=style, linewidth=3, c=color) 
    plt.fill(x, y, alpha=0.1, c=color)


def get_input_ve_model(input_ve_filepath):

    ve = np.genfromtxt(input_ve_filepath, delimiter=',')[2:]
    # "Rcell.exocytosis 0" "G.exocytosis 30" "kt.exocytosis 24" "Npatch.exocytosis 9" 
    # "Nvesicle.exocytosis 21" "Dvesicle.exocytosis 6" "Ninsulin.exocytosis 3" "S.exocytosis 12"
    input_ve_mean = np.array([ve[:,0], ve[:,30], ve[:,24], ve[:,9], ve[:,21], ve[:,6], ve[:,3], ve[:,12]]).transpose()
    input_ve_std = np.array([ve[:,2], ve[:,32], ve[:,26], ve[:,11], ve[:,23], ve[:,8], ve[:,5], ve[:,14]]).transpose()

    return input_ve_mean, input_ve_std


def get_input_pr_model(input_pr_filepath):
    
    pr = np.genfromtxt(input_pr_filepath, delimiter=',')[2:]
    # 'DG.postprandial 18' 'DGd.postprandial 27' 'G.postprandial 15' 'Gb.postprandial 12 '
    # 'S.postprandial 30' 'Sb.postprandial 3' 'I.postprandial 0' 'Y.postprandial 9'
    input_pr_mean = np.array([pr[:,18], pr[:,27], pr[:,15], pr[:,12], pr[:,30], pr[:,3], pr[:,0], pr[:,9]]).transpose()
    input_pr_std = np.array([pr[:,20], pr[:,29], pr[:,17], pr[:,14], pr[:,32], pr[:,5], pr[:,2], pr[:,11]]).transpose()

    return input_pr_mean, input_pr_std


def get_input_pa_model(input_pa_filepath):

    pa = np.genfromtxt(input_pa_filepath, delimiter=',')[2:]
    # 'Spa.pancreas 12' 'Scell.pancreas 3' 'Sis.pancreas 9'
    input_pa_mean = np.array([pa[:,12], pa[:,3], pa[:,9]]).transpose()
    input_pa_std = np.array([pa[:,14], pa[:,5], pa[:,11]]).transpose()

    return input_pa_mean, input_pa_std


def get_input_model(input_ve_filepath,input_pr_filepath,input_pa_filepath):

    input_ve_mean, input_ve_std = get_input_ve_model(input_ve_filepath)
    input_pr_mean, input_pr_std = get_input_pr_model(input_pr_filepath)
    input_pa_mean, input_pa_std = get_input_pa_model(input_pa_filepath)

    return input_ve_mean, input_ve_std, input_pr_mean, input_pr_std, input_pa_mean, input_pa_std



def get_surrogate_ve_model(surrogate_ve_filepath):

    ve = np.genfromtxt(surrogate_ve_filepath, delimiter=',')[2:]
    # "Rcell.exocytosis 0" "G.exocytosis 30" "kt.exocytosis 24" "Npatch.exocytosis 9" 
    # "Nvesicle.exocytosis 21" "Dvesicle.exocytosis 6" "Ninsulin.exocytosis 3" "S.exocytosis 12"
    surrogate_ve_mean = np.array([ve[:,0], ve[:,30], ve[:,24], ve[:,9], ve[:,21], ve[:,6], ve[:,3], ve[:,12]]).transpose()
    surrogate_ve_std = np.array([ve[:,2], ve[:,32], ve[:,26], ve[:,11], ve[:,23], ve[:,8], ve[:,5], ve[:,14]]).transpose()
    
    return surrogate_ve_mean, surrogate_ve_std


def get_surrogate_pr_model(surrogate_pr_filepath):

    pr = np.genfromtxt(surrogate_pr_filepath, delimiter=',')[2:]
    # 'DG.postprandial 18' 'DGd.postprandial 27' 'G.postprandial 15' 'Gb.postprandial 12 '
    # 'S.postprandial 30' 'Sb.postprandial 3' 'I.postprandial 0' 'Y.postprandial 9'
    surrogate_pr_mean = np.array([pr[:,18], pr[:,27], pr[:,15], pr[:,12], pr[:,30], pr[:,3], pr[:,0], pr[:,9]]).transpose()
    surrogate_pr_std = np.array([pr[:,20], pr[:,29], pr[:,17], pr[:,14], pr[:,32], pr[:,5], pr[:,2], pr[:,11]]).transpose()

    return surrogate_pr_mean, surrogate_pr_std


def get_surrogate_pa_model(surrogate_pa_filepath):
    
    pa = np.genfromtxt(surrogate_pa_filepath, delimiter=',')[2:]
    # 'Spa.pancreas 12' 'Scell.pancreas 3' 'Sis.pancreas 9'
    surrogate_pa_mean = np.array([pa[:,12], pa[:,3], pa[:,9]]).transpose()
    surrogate_pa_std = np.array([pa[:,14], pa[:,5], pa[:,11]]).transpose()

    return surrogate_pa_mean, surrogate_pa_std


def get_surrogate_model(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath):
    
    surrogate_ve_mean, surrogate_ve_std = get_surrogate_ve_model(surrogate_ve_filepath)
    surrogate_pr_mean, surrogate_pr_std = get_surrogate_pr_model(surrogate_pr_filepath)
    surrogate_pa_mean, surrogate_pa_std = get_surrogate_pa_model(surrogate_pa_filepath)

    return surrogate_ve_mean, surrogate_ve_std, surrogate_pr_mean, surrogate_pr_std, surrogate_pa_mean, surrogate_pa_std


def get_bnet_metamodel(bnet_filepath):
    
    bnet_result = np.genfromtxt(bnet_filepath, delimiter=',')[2:]

    bnet_ve_mean = np.array([bnet_result[:,45], bnet_result[:,48], bnet_result[:,51], bnet_result[:,54], 
                            bnet_result[:,57], bnet_result[:,60], bnet_result[:,63], bnet_result[:,66]]).transpose()
    bnet_ve_std = np.array([bnet_result[:,47], bnet_result[:,50], bnet_result[:,53], bnet_result[:,56], 
                            bnet_result[:,59], bnet_result[:,62], bnet_result[:,65], bnet_result[:,68]]).transpose()

    bnet_pr_mean = np.array([bnet_result[:,12], bnet_result[:,15], bnet_result[:,18], bnet_result[:,21], 
                            bnet_result[:,24], bnet_result[:,27], bnet_result[:,30], bnet_result[:,33],]).transpose()
    bnet_pr_std = np.array([bnet_result[:,14], bnet_result[:,17], bnet_result[:,20], bnet_result[:,23], 
                            bnet_result[:,26], bnet_result[:,29], bnet_result[:,32], bnet_result[:,35]]).transpose()

    bnet_pa_mean = np.array([bnet_result[:,36], bnet_result[:,39], bnet_result[:,42]]).transpose()
    bnet_pa_std = np.array([bnet_result[:,38], bnet_result[:,41], bnet_result[:,44]]).transpose()

    return  bnet_ve_mean, bnet_ve_std, bnet_pr_mean, bnet_pr_std, bnet_pa_mean, bnet_pa_std


def get_metamodel(meta_filepath):

    meta = np.genfromtxt(meta_filepath, delimiter=',')[2:]

    # "Rcell.exocytosis 33" "G.exocytosis 36" "kt.exocytosis 39" "Npatch.exocytosis 42" 
    # "Nvesicle.exocytosis 45" "Dvesicle.exocytosis 48" "Ninsulin.exocytosis 51" "S.exocytosis 54"
    meta_ve_mean = np.array([meta[:,33], meta[:,36], meta[:,39], meta[:,42], meta[:,45], meta[:,48], meta[:,51], meta[:,54]]).transpose()
    meta_ve_std = np.array([meta[:,35], meta[:,38], meta[:,41], meta[:,44], meta[:,47], meta[:,50], meta[:,53], meta[:,56]]).transpose()

    # 'DG.postprandial 0' 'DGd.postprandial 3' 'G.postprandial 6' 'Gb.postprandial 9 '
    # 'S.postprandial 12' 'Sb.postprandial 15' 'I.postprandial 18' 'Y.postprandial 21'
    meta_pr_mean = np.array([meta[:,0], meta[:,3], meta[:,6], meta[:,9], meta[:,12], meta[:,15], meta[:,18], meta[:,21]]).transpose()
    meta_pr_std = np.array([meta[:,2], meta[:,5], meta[:,8], meta[:,11], meta[:,14], meta[:,17], meta[:,20], meta[:,23]]).transpose()

    # 'Spa.pancreas 24' 'Scell.pancreas 27' 'Sis.pancreas 30'
    meta_pa_mean = np.array([meta[:,24], meta[:,27], meta[:,30]]).transpose()
    meta_pa_std = np.array([meta[:,26], meta[:,29], meta[:,32]]).transpose()

    return meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std


def get_nonopt_metamodel(meta_filepath):

    meta = np.genfromtxt(meta_filepath, delimiter=',')[2:]

    # "Rcell.exocytosis 46" "G.exocytosis 49" "kt.exocytosis 52" "Npatch.exocytosis 55" 
    # "Nvesicle.exocytosis 58" "Dvesicle.exocytosis 61" "Ninsulin.exocytosis 64" "S.exocytosis 67"
    meta_ve_mean = np.array([meta[:,45], meta[:,48], meta[:,51], meta[:,54], meta[:,57], meta[:,60], meta[:,63], meta[:,66]]).transpose()
    meta_ve_std = np.array([meta[:,47], meta[:,50], meta[:,53], meta[:,56], meta[:,59], meta[:,62], meta[:,65], meta[:,68]]).transpose()

    # 'DG.postprandial 10' 'DGd.postprandial 16' 'G.postprandial 19' 'Gb.postprandial 22 '
    # 'S.postprandial 25' 'Sb.postprandial 28' 'I.postprandial 31' 'Y.postprandial 34'
    meta_pr_mean = np.array([meta[:,9], meta[:,15], meta[:,18], meta[:,21], meta[:,24], meta[:,27], meta[:,30], meta[:,33]]).transpose()
    meta_pr_std = np.array([meta[:,11], meta[:,17], meta[:,20], meta[:,23], meta[:,26], meta[:,29], meta[:,32], meta[:,35]]).transpose()

    # 'Spa.pancreas 37' 'Scell.pancreas 40' 'Sis.pancreas 43'
    meta_pa_mean = np.array([meta[:,36], meta[:,39], meta[:,42]]).transpose()
    meta_pa_std = np.array([meta[:,38], meta[:,41], meta[:,44]]).transpose()


    return meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std


def get_nonopt_metamodel_entropy(meta_filepath):

    meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_nonopt_metamodel(meta_filepath)

    meta_pr_entropy = compute_model_entropy(meta_pr_mean, meta_pr_std)
    meta_pa_entropy = compute_model_entropy(meta_pa_mean, meta_pa_std)
    meta_ve_entropy = compute_model_entropy(meta_ve_mean, meta_ve_std)

    meta_entropy = meta_pr_entropy + meta_pa_entropy + meta_ve_entropy

    print('updated surrogate model entropy: ', np.mean(meta_entropy))

    return meta_ve_entropy, meta_pr_entropy, meta_pa_entropy, meta_entropy


def get_surrogate_entropy(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath):

    s_ve_mean, s_ve_std, s_pr_mean, s_pr_std, s_pa_mean, s_pa_std = get_surrogate_model(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)

    s_ve_entropy = compute_model_entropy(s_ve_mean, s_ve_std)
    s_pr_entropy = compute_model_entropy(s_pr_mean, s_pr_std)
    s_pa_entropy = compute_model_entropy(s_pa_mean, s_pa_std)

    surrogate_entropy = s_ve_entropy + s_pr_entropy + s_pa_entropy
    print('surrogate model entropy: ', np.mean(surrogate_entropy))

    return s_ve_entropy, s_pr_entropy, s_pa_entropy, surrogate_entropy



def get_metamodel_entropy(meta_filepath):

    meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_metamodel(meta_filepath)

    meta_pr_entropy = compute_model_entropy(meta_pr_mean, meta_pr_std)
    meta_pa_entropy = compute_model_entropy(meta_pa_mean, meta_pa_std)
    meta_ve_entropy = compute_model_entropy(meta_ve_mean, meta_ve_std)

    meta_entropy = meta_pr_entropy + meta_pa_entropy + meta_ve_entropy

    print('updated surrogate model entropy: ', np.mean(meta_entropy))

    return meta_ve_entropy, meta_pr_entropy, meta_pa_entropy, meta_entropy



def get_model_consistency(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,meta_filepath):

    s_ve_mean, s_ve_std, s_pr_mean, s_pr_std, s_pa_mean, s_pa_std = get_surrogate_model(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)
    meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_metamodel(meta_filepath)

    meta_consistency_ts = [] # consistency over 420 steps

    for timestep in range(len(meta_ve_mean)): 

        meta_consistency_  = []

        for i in range(len(meta_ve_mean[0])):
            meta_ve_con = compute_overlap_of_normal_dist(meta_ve_mean[timestep,i], s_ve_mean[timestep,i], 
                                            meta_ve_std[timestep,i], s_ve_std[timestep,i])

            meta_consistency_ += [meta_ve_con]

        for i in range(len(meta_pr_mean[0])):
            meta_pr_con = compute_overlap_of_normal_dist(meta_pr_mean[timestep,i], s_pr_mean[timestep,i], 
                                            meta_pr_std[timestep,i], s_pr_std[timestep,i])

            meta_consistency_ += [meta_pr_con]

        for i in range(len(meta_pa_mean[0])):
            meta_pa_con = compute_overlap_of_normal_dist(meta_pa_mean[timestep,i], s_pa_mean[timestep,i], 
                                            meta_pa_std[timestep,i], s_pa_std[timestep,i])

            meta_consistency_ += [meta_pa_con]

        meta_consistency_ = np.array(meta_consistency_)
        meta_consistency_ = meta_consistency_[~np.isnan(meta_consistency_)]
        meta_consistency_ts += [np.mean(meta_consistency_)]

    print('model consistency: ', np.mean(meta_consistency_ts))

    return meta_consistency_ts



def get_nonopt_model_consistency(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,meta_filepath):

    s_ve_mean, s_ve_std, s_pr_mean, s_pr_std, s_pa_mean, s_pa_std = get_surrogate_model(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)
    meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_nonopt_metamodel(meta_filepath)

    meta_consistency_ts = [] # consistency over 420 steps

    for timestep in range(len(meta_ve_mean)): 

        meta_consistency_  = []

        for i in range(len(meta_ve_mean[0])):
            meta_ve_con = compute_overlap_of_normal_dist(meta_ve_mean[timestep,i], s_ve_mean[timestep,i], 
                                            meta_ve_std[timestep,i], s_ve_std[timestep,i])

            meta_consistency_ += [meta_ve_con]

        for i in range(len(meta_pr_mean[0])):
            meta_pr_con = compute_overlap_of_normal_dist(meta_pr_mean[timestep,i], s_pr_mean[timestep,i], 
                                            meta_pr_std[timestep,i], s_pr_std[timestep,i])

            meta_consistency_ += [meta_pr_con]

        for i in range(len(meta_pa_mean[0])):
            meta_pa_con = compute_overlap_of_normal_dist(meta_pa_mean[timestep,i], s_pa_mean[timestep,i], 
                                            meta_pa_std[timestep,i], s_pa_std[timestep,i])

            meta_consistency_ += [meta_pa_con]

        meta_consistency_ = np.array(meta_consistency_)
        meta_consistency_ = meta_consistency_[~np.isnan(meta_consistency_)]
        meta_consistency_ts += [np.mean(meta_consistency_)]

    print('model consistency: ', np.mean(meta_consistency_ts))

    return meta_consistency_ts


def get_model_consistenc_pvalue(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath,meta_filepath):

    s_ve_mean, s_ve_std, s_pr_mean, s_pr_std, s_pa_mean, s_pa_std = get_surrogate_model(surrogate_ve_filepath,surrogate_pr_filepath,surrogate_pa_filepath)
    meta_ve_mean, meta_ve_std, meta_pr_mean, meta_pr_std, meta_pa_mean, meta_pa_std = get_metamodel(meta_filepath)

    meta_consistency_ts = [] # consistency over 420 steps

    for timestep in range(len(meta_ve_mean)): 

        meta_consistency_  = []

        for i in range(len(meta_ve_mean[0])):
            meta_ve_con = get_KS_p_value(meta_ve_mean[timestep,i], s_ve_mean[timestep,i], 
                                            meta_ve_std[timestep,i], s_ve_std[timestep,i])

            meta_consistency_ += [meta_ve_con]

        for i in range(len(meta_pr_mean[0])):
            meta_pr_con = get_KS_p_value(meta_pr_mean[timestep,i], s_pr_mean[timestep,i], 
                                            meta_pr_std[timestep,i], s_pr_std[timestep,i])

            meta_consistency_ += [meta_pr_con]

        for i in range(len(meta_pa_mean[0])):
            meta_pa_con = get_KS_p_value(meta_pa_mean[timestep,i], s_pa_mean[timestep,i], 
                                            meta_pa_std[timestep,i], s_pa_std[timestep,i])

            meta_consistency_ += [meta_pa_con]

        meta_consistency_ = np.array(meta_consistency_)
        meta_consistency_ = meta_consistency_[~np.isnan(meta_consistency_)]
        meta_consistency_ts += [np.mean(meta_consistency_)]

    print('model consistency: ', np.mean(meta_consistency_ts))

    return meta_consistency_ts

