from Metamodel import *
# Import libraries:{{{
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributions as dist
import matplotlib.pyplot as plt
from matplotlib import patches
from tqdm import trange
from tqdm import tqdm
from typing import Optional, Union, Callable, Tuple, List, Dict
from torchdiffeq import odeint
from scipy.interpolate import interp1d
from torch.optim.lr_scheduler import ReduceLROnPlateau
from scipy.stats import pearsonr
from matplotlib.pyplot import MultipleLocator
import matplotlib.ticker
from itertools import combinations
from sklearn.feature_selection import mutual_info_regression
from scipy.stats import norm
from torch.distributions import MultivariateNormal
import torch.nn.functional as F
import json
import Models
# }}}

def main():
    outdir = "__glucose_toy_BMM_opt"
    os.makedirs(outdir, exist_ok=True)

    print("Created:", os.path.abspath(outdir))
    print("Exists?", os.path.exists(outdir))


    ###########################################################################
    # NOTE: Do you want to search for the best couplers? 1 == True; 0 == False
    # It will automatically select the best one and use it in the metamodel.
    ###########################################################################
    find_optimal_connecting_var = 0

    f = open("InputModels.json","r")
    data = json.load(f)
    

    mean_scale = 1.0
    obs_noise_scale = 0.15
    method = "MultiScale"

    dt = 0.01
    simulation_time = 8
    steps = int(simulation_time / dt)

    ###########################################################################
    # Trial code block to simplify the main file                              #
    ###########################################################################
    models_data = []
    for model in ["PXR","MM","ADMETAI"]:
        trans_cov_scale =  data[model]["noise_terms"]["trans_cov_scale"]
        emis_cov_scale =  data[model]["noise_terms"]["emis_cov_scale"]
        obs_noise_scale =  data[model]["noise_terms"]["obs_noise_scale"]
        noise_model_type =  data[model]["noise_terms"]["noise_model_type"]

        model_parameters = data[model]["model_parameters"]
        model_dt = data[model]["misc"]["model_dt"]
        steps = int(simulation_time / model_dt)
        model_function = getattr(Models, data[model]["misc"]["model_function_name"])
        # simulate input model ODE
        df_model, plot_model, model_data_tensor = simulate_model(model_function, model_parameters, steps, simulation_time)

        model_std = (model_data_tensor * 0.15)
        ylabels = [f"{s} [{model_parameters['units'][i]}]" for i, s in enumerate(model_parameters["species"])]
        fig = plot_inputmodel(model_data_tensor, model_std, dt=model_dt, sim_time=simulation_time, name=ylabels)
        figname = os.path.join(outdir, data[model]["misc"]["input_model_plot"])
        fig.tight_layout()
        fig.savefig(figname)

        
        obs_model, t_model, df_model = get_observations(
            rate_eqs=model_function, parameters=model_parameters, dt=model_dt,
            sim_time=simulation_time, noise_scale=obs_noise_scale, kwargs=None,
        )
        obs_model[:, :, 0] *= mean_scale
        fig = plot_inputmodel(model_data_tensor, model_std, dt=model_dt, sim_time=simulation_time, name=ylabels, observables=obs_model)
        figname = os.path.join(outdir, data[model]["misc"]["input_model_with_obs_plot"])
        fig.tight_layout()
        fig.savefig(figname)

        # surrogate, surrogate_states, surrogate_df, param_names = generate_surrogate(
        #     observables=obs_model, model_parameters=model_parameters,
        #     ode_fn=model_function, sim_time=simulation_time, dt=model_dt,
        #     noise_model_type='time-variant'
        # )

        surrogate, _, _, param_names = generate_surrogate(
            observables=obs_model, model_parameters=model_parameters,
            ode_fn=model_function, sim_time=simulation_time, dt=model_dt,
            noise_model_type='time-variant'
        )
        # print(type(surrogate_states))
        # input()
        # fig = plot_inputmodel(
        #     input_result=surrogate_states,
        #     input_std=obs_model[:, :, 1],
        #     dt=model_dt,
        #     sim_time=simulation_time,
        #     name=ylabels,
        #     observables=obs_model
        # )
        # figname = os.path.join(outdir, data[model]["misc"]["surrogate_with_observables_plot"])
        # fig.tight_layout()
        # fig.savefig(figname)

        models_data_current = \
        {
            'obs_model': obs_model,
            'surrogate': surrogate,
            'params': model_parameters,
            'ode_fn': model_function,
            'dt': model_dt,
            'transition_cov_scale': trans_cov_scale,
            'emission_cov_scale': emis_cov_scale,
            'observation_noise_scale': obs_noise_scale,
            'noise_model_type': noise_model_type,
            'param_names': param_names
        }
        models_data.append(models_data_current)

    ##########################################################################

    if find_optimal_connecting_var:
        connecting_vars = None
    else:
        connecting_vars = [((0, 1),('$C^{PXR}_{CYP}$','$C^{MM}_{CYP}$'),('physical')),
                           ((1, 0),('$C^{MM}_{SUB}$','$C^{PXR}_{RIF}$'),('physical')),
                           ((1, 2),('$C^{MM}_{SUB}$','$C^{ADMETAI}_{SUB}$'),('statistical'))]
        print("connecting_vars =", connecting_vars)

    ###########################################################################
    # NOTE: This is where the action is... Here, we build the metamodel and
    # specify if we want to optimize the coupling weights.
    ###########################################################################
    bayes_mm, coupled_states, obs_dfs, coupled_dfs, figs, connecting_vars,\
            processed_obs, covariances_over_time = build_metamodel(
        models_data=models_data,
        sim_time=simulation_time,
        connecting_vars=connecting_vars
    )



    obs_model_human_resampled = processed_obs[0]
    obs_model_rabbit_resampled = processed_obs[1]

    #print("Optimized Coupling Weights:", optimized_weights)
#    print("optimized_params = ", optimized_params)
    states_a, states_b, states_c = coupled_states
    data_df_a, data_df_b, data_df_c = obs_dfs
    df_a, df_b, df_c = coupled_dfs
    times = df_a["time"].to_numpy()
    universal_dt = times[1] - times[0]

    ###########################################################################
    # NOTE: Obviously, you will need to modify the plotting function for a
    # different example. And you might not have a "ground-truth".
    ###########################################################################
    fig = plot_groundtruth_with_metamodel(
        metamodel_a=states_a,
        metamodel_b=states_b,
        metamodel_std_a=np.stack([data_df_a[f"{name}_std"] for name in data["PXR"]["model_parameters"]["species"]], axis=1),
        metamodel_std_b=np.stack([data_df_b[f"{name}_std"] for name in data["MM"]["model_parameters"]["species"]], axis=1),
        model_a_params=data["PXR"]["model_parameters"],
        model_b_params=data["MM"]["model_parameters"],
        model_PXR_dt=data["PXR"]["misc"]["model_dt"],
        model_human_dt=data["PXR"]["misc"]["model_dt"],
        model_rabbit_dt=data["PXR"]["misc"]["model_dt"],
        sim_time=simulation_time,
        observables_a=obs_model_human_resampled,
        observables_b=obs_model_rabbit_resampled
    )
    figname = os.path.join(outdir, "groundtruth_with_metamodel.png")
    fig.tight_layout()
    fig.savefig(figname)
    
    compute_clearance_from_metamodel(
        metamodel_a=states_b,
        metamodel_std_a=np.stack([data_df_b[f"{name}_std"] for name in data["MM"]["model_parameters"]["species"]], axis=1),
        model_a_params=data["MM"]["model_parameters"],
        model_dt=data["MM"]["misc"]["model_dt"],
        sim_time=simulation_time)
    
    
if __name__ == "__main__":
    main()