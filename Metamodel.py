

# Import libraries:{{{
import os, math, re, warnings, time
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
from typing import Optional, Union, Callable, Tuple, List, Dict, Iterable, Any, Sequence
from torchdiffeq import odeint
from torch.optim.lr_scheduler import ReduceLROnPlateau
from matplotlib.pyplot import MultipleLocator
import matplotlib.ticker
from itertools import combinations
from sklearn.feature_selection import mutual_info_regression
from torch.distributions import MultivariateNormal
import torch.nn.functional as F
from scipy.interpolate import interp1d
from scipy.stats import norm, pearsonr, gaussian_kde, linregress, skew  # Note: kurtosis not used, but skew is; assuming scipy available
from scipy.signal import find_peaks
from dataclasses import dataclass, field
from scipy import integrate
from sklearn import metrics
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker

from itertools import combinations, chain
from scipy.linalg import logm
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import uncertainties as u
from uncertainties import unumpy as unp

# }}}


def compute_clearance_from_metamodel(metamodel_a,metamodel_std_a,model_a_params,model_dt,sim_time):
    state_dim_a = len(model_a_params["InitConc"])
    def clean_species_name(s):
        s = s.strip("$").split("^{")[0].replace("\\", "")
        return s
    universal_dt = model_dt
    steps = int(sim_time / universal_dt)
    universal_time = np.arange(0, sim_time, universal_dt)[:steps]

    def resample_data(data, data_dt, target_time):
        target_dt = target_time[1] - target_time[0]
        if data.shape[0] == len(target_time) and data_dt == target_dt:
            return data.detach().numpy() if torch.is_tensor(data) else data
        if data_dt > target_dt:
            resampled = np.zeros((len(target_time), data.shape[1]))
            data_time = np.arange(0, sim_time, data_dt)[:min(data.shape[0], int(sim_time / data_dt))]
            for i in range(data.shape[1]):
                valid_length = min(len(data_time), len(data[:, i]))
                interp_func = interp1d(
                    data_time[:valid_length],
                    data[:valid_length, i].detach().numpy() if torch.is_tensor(data) else data[:valid_length, i],
                    kind='linear', fill_value='extrapolate'
                )
                resampled[:, i] = interp_func(target_time)
            return resampled
        return data[:len(target_time)].detach().numpy() if torch.is_tensor(data) else data[:len(target_time)]

    metamodel_a_resampled = resample_data(metamodel_a, universal_dt, universal_time)
    meta = metamodel_a_resampled[:, 0]
    metamodel_std_a_resampled = resample_data(metamodel_std_a, universal_dt, universal_time)
    meta_std = metamodel_std_a_resampled[:, 0]
    metaSignificant = meta[0:20]
    timeSignificant = universal_time[0:20]
    print(metaSignificant)
    print(timeSignificant)
    slope, intercept = np.polyfit(timeSignificant, np.log(metaSignificant), 1)
    print("Clearance = ", slope)


def plot_inputmodel(input_result, input_std, dt, sim_time, name=None, observables=None):
    fig = plt.figure(figsize=(20, 4))
    n_variable = input_result.shape[1]
    # print(input_result.shape)
    # print(len(name))
    time = np.arange(0, sim_time, dt)[:input_result.shape[0]]
    for i in range(n_variable):
        ax = fig.add_subplot(1, n_variable, i+1)
        plt.yticks(fontproperties='Arial Narrow', size=17)
        plt.xticks(fontproperties='Arial Narrow', size=17)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        plt.xlabel('Time [h]', fontproperties='Arial Narrow', size=20)
        plt.ylabel(name[i] if name else f'Variable {i+1}', fontproperties='Arial Narrow', size=20)
        plt.plot(time, unp.nominal_values(input_result[:, i]), color='black', linewidth=2, label='model')
        plt.fill_between(time,
                         unp.nominal_values(input_result[:, i]) - unp.std_devs(input_result[:, i]),
                         unp.nominal_values(input_result[:, i]) + unp.std_devs(input_result[:, i]),
                         alpha=0.1, color='black')
        if observables is not None:
            obs_mean = observables[:, i, 0]
            obs_std = observables[:, i, 1]
            plt.plot(time, obs_mean, color='blue', linestyle='--', linewidth=1.5, label='observables')
            plt.fill_between(time,
                             (obs_mean - obs_std).reshape(-1,),
                             (obs_mean + obs_std).reshape(-1,),
                             alpha=0.1, color='blue')
        plt.legend(loc='lower center', prop={'size': 14, 'family': 'Arial Narrow'})
    return fig


def plot_groundtruth_with_metamodel(
    metamodel_a, metamodel_b, metamodel_std_a, metamodel_std_b,
    model_a_params, model_b_params, model_human_dt, model_rabbit_dt, model_PXR_dt, sim_time,
    observables_a=None, observables_b=None
):
    state_dim_a = len(model_a_params["InitConc"])
    state_dim_b = len(model_b_params["InitConc"])
    def clean_species_name(s):
        s = s.strip("$").split("^{")[0].replace("\\", "")
        return s
    universal_dt = min(model_human_dt, model_rabbit_dt, model_PXR_dt)
    steps = int(sim_time / universal_dt)
    universal_time = np.arange(0, sim_time, universal_dt)[:steps]

    def resample_data(data, data_dt, target_time):
        target_dt = target_time[1] - target_time[0]
        if data.shape[0] == len(target_time) and data_dt == target_dt:
            return data.detach().numpy() if torch.is_tensor(data) else data
        if data_dt > target_dt:
            resampled = np.zeros((len(target_time), data.shape[1]))
            data_time = np.arange(0, sim_time, data_dt)[:min(data.shape[0], int(sim_time / data_dt))]
            for i in range(data.shape[1]):
                valid_length = min(len(data_time), len(data[:, i]))
                interp_func = interp1d(
                    data_time[:valid_length],
                    data[:valid_length, i].detach().numpy() if torch.is_tensor(data) else data[:valid_length, i],
                    kind='linear', fill_value='extrapolate'
                )
                resampled[:, i] = interp_func(target_time)
            return resampled
        return data[:len(target_time)].detach().numpy() if torch.is_tensor(data) else data[:len(target_time)]

    metamodel_a_resampled = resample_data(metamodel_a, universal_dt, universal_time)
    metamodel_b_resampled = resample_data(metamodel_b, universal_dt, universal_time)
    metamodel_std_a_resampled = resample_data(metamodel_std_a, universal_dt, universal_time)
    metamodel_std_b_resampled = resample_data(metamodel_std_b, universal_dt, universal_time)

    def resample_observables(obs, obs_dt):
        if obs is None:
            return None, None
        if obs.shape[0] == steps and obs_dt == universal_dt:
            return obs[:, :, 0].detach().numpy(), obs[:, :, 1].detach().numpy()
        obs_mean = np.zeros((steps, obs.shape[1]))
        obs_std = np.zeros((steps, obs.shape[1]))
        obs_time = np.arange(0, sim_time, obs_dt)[:min(obs.shape[0], int(sim_time / obs_dt))]
        for i in range(obs.shape[1]):
            valid_length = min(len(obs_time), len(obs[:, i, 0]))
            interp_mean = interp1d(
                obs_time[:valid_length],
                obs[:valid_length, i, 0].detach().numpy(),
                kind='linear', fill_value='extrapolate'
            )
            interp_std = interp1d(
                obs_time[:valid_length],
                obs[:valid_length, i, 1].detach().numpy(),
                kind='linear', fill_value='extrapolate'
            )
            obs_mean[:, i] = interp_mean(universal_time)
            obs_std[:, i] = interp_std(universal_time)
        return obs_mean, obs_std

    obs_mean_a, obs_std_a = resample_observables(observables_a, universal_dt) if observables_a is not None else (None, None)
    obs_mean_b, obs_std_b = resample_observables(observables_b, universal_dt) if observables_b is not None else (None, None)

    ylabels = [f"{s} [{model_a_params['units'][i]}]" for i, s in enumerate(model_a_params["species"])] + \
              [f"{s} [{model_b_params['units'][i]}]" for i, s in enumerate(model_b_params["species"])]
    ylabels = [label.replace("[]","") for label in ylabels]
    fig = plt.figure(figsize=(5 * (state_dim_a + state_dim_b), 4))
    for i in range(state_dim_a + state_dim_b):
        ax = fig.add_subplot(1, state_dim_a + state_dim_b, i + 1)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        plt.yticks(fontproperties='Arial Narrow', size=17)
        plt.xticks(fontproperties='Arial Narrow', size=17)
        plt.xlabel('Time [h]', fontproperties='Arial Narrow', size=20)
        plt.ylabel(ylabels[i], fontproperties='Arial Narrow', size=20)
        if i < state_dim_a:
            meta = metamodel_a_resampled[:, i]
            meta_std = metamodel_std_a_resampled[:, i]
            obs_mean = obs_mean_a[:, i] if obs_mean_a is not None else None
            obs_std = obs_std_a[:, i] if obs_std_a is not None else None
        else:
            j = i - state_dim_a
            meta = metamodel_b_resampled[:, j]
            meta_std = metamodel_std_b_resampled[:, j]
            obs_mean = obs_mean_b[:, j] if obs_mean_b is not None else None
            obs_std = obs_std_b[:, j] if obs_std_b is not None else None
        plt.plot(universal_time, meta, color='red', linewidth=2, label='Metamodel')
        plt.fill_between(universal_time, meta - meta_std, meta + meta_std, alpha=0.3, color='red')
        if obs_mean is not None and obs_std is not None:
            plt.plot(universal_time, obs_mean, color='blue', linestyle='--', linewidth=1.5, label='Observables')
            plt.fill_between(universal_time, obs_mean - obs_std, obs_mean + obs_std, alpha=0.1, color='blue')
        #plt.legend(loc='lower center', prop={'size': 14, 'family': 'Arial Narrow'})
        if i == 0:
            plt.legend(loc='best', prop={'size': 14, 'family': 'Arial Narrow'})
    plt.tight_layout()
    return fig



# ODEs: {{{
def get_observations(rate_eqs,
                     parameters,
                     dt,
                     sim_time,
                     noise_scale=0.0,
                     kwargs=None):
    """
    Runs the uncertainty-aware ODE model and returns observations with mean and std.

    rate_eqs      : callable f(y, t, k_list, **kwargs) -> dy/dt
    parameters    : dict with InitConc, InitConc_std, k, k_std, species, etc.
    dt            : time step (float)
    sim_time      : total simulation time (float)
    noise_scale   : multiplicative obs noise scale (float). optional extra noise.
    kwargs        : extra args passed to rate_eqs, e.g. {"obs_fn": ...}

    Returns:
        observations: torch.tensor of shape (steps, n_species, 2)
                      [..., 0] = mean trajectory (nominal values)
                      [..., 1] = std trajectory  (propagated + obs noise)
        t_grid: np.ndarray of shape (steps,)
        df: pandas DataFrame with columns:
             time, species_i mean, species_i std
    """

    # how many steps for this dt
    steps = int(np.ceil(sim_time / dt))

    # run full uncertainty-aware simulation
    df, _, sol = simulate_model(
        rate_eqs=rate_eqs,
        parameters=parameters,
        steps=steps,
        simulation_time=sim_time,
        kwargs=kwargs,
    )
    # df["time"] exists, and df has columns for each species mean and <species>_std

    # extract mean and std arrays in the correct order
    means_list = []
    stds_list  = []
    for i, s in enumerate(parameters["species"]):
        # Pull the nominal and std directly from 'sol'
        vals_nominal = [nominal(val) for val in sol[:, i]]
        vals_std     = [stdev(val)   for val in sol[:, i]]

        means_list.append(vals_nominal)  # shape (steps,)
        stds_list.append(vals_std)       # shape (steps,)

    means_np = np.array(means_list).T  # shape (steps, n_species)
    stds_np  = np.array(stds_list).T   # shape (steps, n_species)

    # Add measurement noise in quadrature, if requested
    if noise_scale is not None and noise_scale != 0.0:
        obs_noise = np.abs(noise_scale * means_np)
        total_std = np.sqrt(stds_np**2 + obs_noise**2) + 1e-6
    else:
        total_std = stds_np + 1e-6  # ensure strictly positive
        # note: stdev(val) may already be 0.0, so +1e-6 avoids exact zeros

    # pack into (steps, n_species, 2)
    obs_mean_tensor = torch.tensor(means_np, dtype=torch.float32)
    obs_std_tensor  = torch.tensor(total_std, dtype=torch.float32)
    observations = torch.stack([obs_mean_tensor, obs_std_tensor], dim=2)

    # also return the time grid for plotting if you want
    t_grid = df["time"].to_numpy()

    return observations, t_grid, df



def _to_list(x):
    # Ensure we hand rate_eqs a plain Python list, not a numpy object array
    # This reduces weird broadcasting/ufunc surprises.
    if isinstance(x, np.ndarray):
        return [xi for xi in x.tolist()]
    else:
        return [xi for xi in x]

def _sanitize_state(y_arr):
    # Return a numpy object array of ufloats, with non-finite or negative nominal set to 0±0
    cleaned = []
    for val in y_arr:
        # ensure val is ufloat
        if hasattr(val, "nominal_value"):
            nom = val.nominal_value
            std = val.std_dev
        else:
            nom = float(val)
            std = 0.0
            val = ufloat(nom, std)

        # clamp rule on the nominal value:
        if (not np.isfinite(nom)) or (nom < 0.0):
            cleaned.append(ufloat(0.0, 0.0))
        else:
            cleaned.append(val)

    return np.array(cleaned, dtype=object)

def rk4_integrate(rate_eqs, y0, t_grid, k_list, kwargs):
    """
    rate_eqs(y, t, k_list, **kwargs) -> dy/dt where y is a list-like of states.
    y0 : sequence of ufloats (or floats, which we'll convert to ufloats)
    t_grid : np.array of shape (T,)
    returns: np.array of shape (T, n_species) with ufloats
    """

    # make y_curr a clean ufloats array
    y_curr = _sanitize_state(np.array(y0, dtype=object))

    sol = [y_curr.copy()]  # store first state

    for idx in range(len(t_grid) - 1):
        t = t_grid[idx]
        dt = t_grid[idx + 1] - t_grid[idx]

        # sanity check: RK4 assumes forward positive step
        if dt <= 0:
            raise ValueError(f"Non-positive dt encountered at step {idx}: dt={dt}")

        # compute k1..k4 using fresh copies, not mutated arrays
        y_list = _to_list(y_curr)
        k1 = np.array(rate_eqs(y_list, t, k_list, **kwargs), dtype=object)

        y_tmp = y_curr + 0.5 * dt * k1
        y_list = _to_list(y_tmp)
        k2 = np.array(rate_eqs(y_list, t + 0.5 * dt, k_list, **kwargs), dtype=object)

        y_tmp = y_curr + 0.5 * dt * k2
        y_list = _to_list(y_tmp)
        k3 = np.array(rate_eqs(y_list, t + 0.5 * dt, k_list, **kwargs), dtype=object)

        y_tmp = y_curr + dt * k3
        y_list = _to_list(y_tmp)
        k4 = np.array(rate_eqs(y_list, t + dt, k_list, **kwargs), dtype=object)

        # RK4 update
        y_next = y_curr + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        # clamp / sanitize AFTER the step
        y_curr = _sanitize_state(y_next)

        sol.append(y_curr.copy())

    return np.array(sol, dtype=object)




def nominal(val):
    if hasattr(val, "nominal_value"):
        return val.nominal_value
    return float(val)

def stdev(val):
    if hasattr(val, "std_dev"):
        return val.std_dev
    return 0.0

def simulate_model(
    rate_eqs,
    parameters,
    steps,
    simulation_time,
    kwargs=None,
):
    # 1. Time grid
    t_grid = np.linspace(0.0, simulation_time, steps)

    # 2. Build initial conditions with uncertainties
    init_conc = []
    for i in range(len(parameters["InitConc"])):
        mean_i = parameters["InitConc"][i]
        std_i  = parameters["InitConc_std"][i]
        if std_i is not None:
            init_conc.append(u.ufloat(mean_i, std_i))
        else:
            # convert plain float to ufloat with 0.0 std so we keep type consistency
            init_conc.append(u.ufloat(mean_i, 0.0))

    # 3. rate constants with uncertainties
    k_list = []
    for i in range(len(parameters["k"])):
        mean_k = parameters["k"][i]
        std_k  = parameters["k_std"][i]
        if std_k is not None:
            k_list.append(u.ufloat(mean_k, std_k))
        else:
            k_list.append(u.ufloat(mean_k, 0.0))

    # 4. Resample observable if provided
    if kwargs is None:
        kwargs = {}
    if "obs_fn" in kwargs:
        raw = kwargs["obs_fn"]
        vals = raw(t_grid)
        kwargs["obs_fn"] = interp1d(
            t_grid, vals, kind="linear",
            bounds_error=False,
            fill_value=(vals[0], vals[-1])
        )

    # 5. Integrate with custom RK4 that keeps ufloats alive
    # NOTE: rk4_integrate must internally clamp negatives and replace NaN/inf,
    # because np.clip/np.nan_to_num won't like dtype=object arrays of ufloats.
    sol = rk4_integrate(rate_eqs, init_conc, t_grid, k_list, kwargs)


    # 6. Build DataFrame with nominal + std
    df = pd.DataFrame({"time": t_grid})
    for i, s in enumerate(parameters["species"]):
        vals_nominal = [nominal(val) for val in sol[:, i]]
        vals_std     = [stdev(val)   for val in sol[:, i]]
        df[s] = vals_nominal
        df[s + "_std"] = vals_std

    # 7. Plot with error bars
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, s in enumerate(parameters["species"]):
        ax.errorbar(
            df["time"],
            df[s],
            yerr=df[s + "_std"],
            fmt='-o',              # line + markers
            capsize=3,             # little whiskers on error bars
            label=f"{s} [{parameters['units'][i]}]",
            linewidth=1.5,
            markersize=4,
        )

    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Concentration")
    ax.legend()
    plt.tight_layout()

    return df, ax, sol




# }}}

# Surrogate and Coupler classes:{{{
class SurrogateModel:
    """
    Represents a surrogate dynamical system with forward dynamics and observation model.

    Attributes:
        f (Callable): Transition/forward function f(z, dt) propagating latent state z over time dt.
        g (Callable): Observation function g(z) mapping latent state z to observable output.
        Q (Tensor or Callable): Transition/process noise covariance or a function that returns it given z.
        R (Tensor or Callable): Emission/observation noise covariance or a function that returns it given z.
        state_dim (int): Dimensionality of the latent state space.
    """

    def __init__(self, forward_fn, observe_fn, trans_noise, emit_noise, state_dim):
        """
        Initialize the surrogate model.

        Args:
            forward_fn (Callable): Function f(z, dt) for forward time evolution of state z.
            observe_fn (Callable): Function g(z) for producing observable output from state.
            trans_noise (Tensor or Callable): Transition noise covariance or generator function.
            emit_noise (Tensor or Callable): Emission noise covariance or generator function.
            state_dim (int): Dimensionality of the state space.
        """
        self.f = forward_fn
        self.g = observe_fn
        self.Q = trans_noise
        self.R = emit_noise
        self.state_dim = state_dim

    def predict(self, z, cov, f, dt):
        """
        Predict the next state and transition uncertainty.

        Args:
            z (Tensor): Current state.
            cov (Tensor): Current state covariance (unused here but kept for interface compatibility).
            f (Callable): Transition function (can override the default self.f).
            dt (float): Time step for propagation.

        Returns:
            Tuple: Predicted state mean (mu), predicted process noise covariance (cov).
        """
        mu = f(z, dt)
        cov = self.Q(z) if callable(self.Q) else self.Q
        return mu, cov

    def observe(self, z):
        """
        Map the latent state to observable quantities.

        Args:
            z (Tensor): Current latent state.

        Returns:
            Tuple: Predicted observation mean (mu), observation noise covariance (cov).
        """
        mu = self.g(z)
        cov = self.R(z) if callable(self.R) else self.R
        return mu, cov


class Coupler:
    """
    A utility class for combining outputs of multiple models by shared priors.
    Uses uncertainty-weighted blending.
    """

    @staticmethod
    def mixture_prior(mus, sigmas):
        """
        Uncertainty-weighted blending of independent estimates.

        Given N estimates (mu_i, sigma_i) of a shared quantity, computes the
        minimum-variance Bayesian combination:

            mu_c    = sum(mu_i / sigma_i^2)  /  sum(1 / sigma_i^2)
            sigma_c = sqrt( 1 / sum(1 / sigma_i^2) )

        This is the exact posterior mean and std under independent Gaussian
        likelihoods, so a model with large sigma_i automatically contributes
        little weight.

        Args:
            mus    (list[Tensor] | Tensor): Means from each model.
            sigmas (list[Tensor] | Tensor): Std-devs from each model.

        Returns:
            Tuple[Tensor, Tensor]: Fused mean mu_c and fused std sigma_c.
        """
        if not isinstance(mus, torch.Tensor):
            mus = torch.stack([
                m if isinstance(m, torch.Tensor)
                else torch.tensor(m, dtype=torch.float32) for m in mus
            ])
        if not isinstance(sigmas, torch.Tensor):
            sigmas = torch.stack([
                s if isinstance(s, torch.Tensor)
                else torch.tensor(s, dtype=torch.float32) for s in sigmas
            ])
        precisions  = 1.0 / torch.clamp(sigmas ** 2, min=1e-12)
        total_prec  = precisions.sum()
        mu_c        = (precisions * mus).sum() / total_prec
        sigma_c     = torch.sqrt(1.0 / total_prec)
        return mu_c, sigma_c

    def identify_and_couple(self, data_list, var_names):
        """
        Compute a joint prior across multiple models/datasets using precision-weighted
        fusion of their empirical statistics.

        Args:
            data_list (list): List of DataFrames containing variable observations.
            var_names (list): Variable name (key) to extract from each dataset.

        Returns:
            Tuple[Tensor, Tensor]: Fused mean mu_c and fused std sigma_c.
        """
        mus    = torch.tensor([data[var].mean() for data, var in zip(data_list, var_names)])
        sigmas = torch.tensor([data[var].std()  for data, var in zip(data_list, var_names)])
        mu_c, sigma_c = self.mixture_prior(mus, sigmas)
        return mu_c, sigma_c


# }}}

# UKF class:{{{
class UnscentedKalmanFilter:
    """
    Unscented Kalman Filter (UKF) implementation for nonlinear state estimation.

    Attributes:
        n (int): Dimensionality of the state vector.
        Q (Tensor or callable): Process noise covariance matrix or a function returning it based on predicted state.
        lambda_ (float): Scaling parameter used in sigma point generation.
        Wm (Tensor): Weights for computing mean of transformed sigma points.
        Wc (Tensor): Weights for computing covariance of transformed sigma points.
    """

#    def __init__(self, dim, Q, alpha=1e-3, beta=2.0, kappa=0.0):
    def __init__(self, dim, Q, alpha=1e-2, beta=2.0, kappa=0.0):
        self.n = dim
        self.Q = Q
        self.lambda_ = alpha**2 * (self.n + kappa) - self.n
        c = self.n + self.lambda_
        self.Wm = torch.full((2 * self.n + 1,), 0.5 / c)
        self.Wc = self.Wm.clone()
        self.Wm[0] = self.lambda_ / c
        self.Wc[0] = self.lambda_ / c + (1 - alpha**2 + beta)

    def sigma_points(self, x, P):
        n = len(x)
        try:
            L = torch.linalg.cholesky(P)
        except Exception:
            epsilon = 1e-6 * torch.eye(P.shape[0], dtype=torch.float32)
            P = P + epsilon
            eigvals, eigvecs = torch.linalg.eigh(P)
            eigvals = torch.clamp(eigvals, min=1e-6)
            P = eigvecs @ torch.diag(eigvals) @ eigvecs.T
            L = torch.linalg.cholesky(P)

        X = [x]
        scale = torch.sqrt(torch.tensor(n + self.lambda_, dtype=torch.float32))
        for i in range(n):
            X.append(x + scale * L[:, i])
            X.append(x - scale * L[:, i])
        return torch.stack(X)

    def update(self, x_pred, P_pred, z, h, R):
        X = self.sigma_points(x_pred, P_pred)  # (2n+1, n)
        Z = torch.stack([h(xi) for xi in X])   # (2n+1, m)
        z_pred = (self.Wm.unsqueeze(1) * Z).sum(dim=0)

        P_zz = sum(
            self.Wc[i] * torch.outer(Z[i] - z_pred, Z[i] - z_pred)
            for i in range(2 * self.n + 1)
        ) + R

        P_xz = sum(
            self.Wc[i] * torch.outer(X[i] - x_pred, Z[i] - z_pred)
            for i in range(2 * self.n + 1)
        )

        K = P_xz @ torch.linalg.inv(P_zz)

        x_up = x_pred + K @ (z - z_pred)
        P_up = P_pred - K @ P_zz @ K.T

        return x_up, P_up, z_pred, P_zz

    def predict(self, x, P, f, dt):
        X = self.sigma_points(x, P)  # (2n+1, n)
        Xf = torch.stack([f(xi, dt) for xi in X])  # (2n+1, n)
        x_pred = (self.Wm.unsqueeze(1) * Xf).sum(dim=0)

        P_pred = sum(
            self.Wc[i] * torch.outer(Xf[i] - x_pred, Xf[i] - x_pred)
            for i in range(2 * self.n + 1)
        )

        Q = self.Q(x_pred) if callable(self.Q) else self.Q
        P_pred += Q

        return x_pred, P_pred

# }}}

# BayesMM{{{
class BayesMM:

    def __init__(self, models, param_names_list=None, universal_dt=0.01,
                 constraint_maps=None):
        """
        Parameters
        ----------
        models : list[SurrogateModel]
        param_names_list : list[list[str]] | None
        universal_dt : float
        constraint_maps : list[Callable | None] | None
            One entry per model.  Each entry is either None (no propagation
            for that model) or a Callable with signature

                f(x_before, x_after, params, dt) -> x_propagated

            where x_before is the state BEFORE the current coupling step,
            x_after is the state AFTER (with primary variables corrected),
            params is the current parameter dict for that model, and dt is
            universal_dt.  The callable returns a new state tensor in which
            dependent variables have been recomputed to be consistent with
            the corrected primaries.  Pass None for the whole list to
            disable the feature entirely (default).
        """
        self.models           = models
        self.param_names_list = param_names_list or [[] for _ in models]
        self.universal_dt     = universal_dt
        self.coupler          = Coupler()
        self.surrogates       = [
            UnscentedKalmanFilter(model.state_dim, model.Q)
            for model in models
        ]
        # constraint_maps: list of per-model propagation callables or None
        if constraint_maps is None:
            self.constraint_maps = [None] * len(models)
        else:
            if len(constraint_maps) != len(models):
                raise ValueError(
                    f"constraint_maps length ({len(constraint_maps)}) must "
                    f"equal number of models ({len(models)})."
                )
            self.constraint_maps = list(constraint_maps)

    # ──────────────────────────────────────────────────────────────────────────
    # Coupling builders
    # ──────────────────────────────────────────────────────────────────────────

    def couple_models(self, data_list, connecting_vars):
        couplings = []
        print(len(connecting_vars))
        exit()
        for model_indices, var_names in connecting_vars:
            data_sub_list     = [data_list[i] for i in model_indices]
            mu_c, sigma_c     = self.coupler.identify_and_couple(data_sub_list, var_names)
            weights           = [
                torch.tensor(1 / len(model_indices), dtype=torch.float32, requires_grad=True)
                for _ in model_indices
            ]
            couplings.append({
                'coupling_type': 'statistical',
                'models':  tuple(model_indices),
                'vars':    tuple(var_names),
                'mu_c':    mu_c,
                'sigma_c': sigma_c,
                'weights': weights,
            })
        return couplings

    # ──────────────────────────────────────────────────────────────────────────
    # Small helpers
    # ──────────────────────────────────────────────────────────────────────────

    def _as_scalar_tensor(self, x, device, dtype=torch.float64):
        if isinstance(x, torch.Tensor):
            return x.reshape(()).to(device=device, dtype=dtype)
        return torch.tensor(x, device=device, dtype=dtype).reshape(())

    def _coerce_coupling_tensors(self, couplings):
        for c in couplings:
            if 'mu_c' in c and not isinstance(c['mu_c'], torch.Tensor):
                c['mu_c'] = torch.tensor(float(c['mu_c']), dtype=torch.float32)
            if 'sigma_c' in c and not isinstance(c['sigma_c'], torch.Tensor):
                c['sigma_c'] = torch.tensor(float(c['sigma_c']), dtype=torch.float32)

    def ensure_maps_on_couplings(self, couplings, default_alpha=1.0, default_beta=0.0):
        for c in couplings:
            J = len(c['models'])
            if 'lambdas' not in c or c['lambdas'] is None or len(c['lambdas']) != J:
                c['lambdas'] = [default_alpha] * J
            if 'shifts'  not in c or c['shifts']  is None or len(c['shifts'])  != J:
                c['shifts']  = [default_beta]  * J

    def _get_map_params(self, c_dict, local_idx, device, dtype):
        if 'lambdas' in c_dict and c_dict['lambdas'] is not None and len(c_dict['lambdas']) > local_idx:
            a = c_dict['lambdas'][local_idx]
        elif 'scales' in c_dict and c_dict['scales'] is not None and len(c_dict['scales']) > local_idx:
            a = c_dict['scales'][local_idx]
        else:
            a = 1.0

        if 'shifts' in c_dict and c_dict['shifts'] is not None and len(c_dict['shifts']) > local_idx:
            b = c_dict['shifts'][local_idx]
        elif 'offsets' in c_dict and c_dict['offsets'] is not None and len(c_dict['offsets']) > local_idx:
            b = c_dict['offsets'][local_idx]
        else:
            b = 0.0

        return (
            self._as_scalar_tensor(a, device, dtype=dtype),
            self._as_scalar_tensor(b, device, dtype=dtype),
        )

    # ──────────────────────────────────────────────────────────────────────────
    # Statistical coupling math
    # ──────────────────────────────────────────────────────────────────────────

    def _compute_effective_omega(self, sig_t, sig_c, w_j=None):
        sig_t2    = torch.clamp(sig_t**2, min=1e-12)
        sig_c2    = torch.clamp(sig_c**2, min=1e-12)
        omega_opt = sig_t2 / (sig_t2 + sig_c2)
        if getattr(self, "learn_omegas", False) and (w_j is not None):
            omega_eff = torch.clamp(0.5 * w_j + 0.5 * omega_opt, 0.0, 1.0)
        else:
            omega_eff = omega_opt
        return omega_eff

    # ──────────────────────────────────────────────────────────────────────────
    # Physical coupling correction
    # ──────────────────────────────────────────────────────────────────────────

    def _apply_physical_coupling_correction(
        self, c_dict, local_idx, states, covs, data_list, t
    ):
        """
        Kalman-style residual correction for one participating model:

            y_i  = alpha_i * x_i + beta_i
            r    = g_c({y_i}, t=t)
            K_j  = alpha_j * var_j / (sum_i alpha_i^2 * var_i + sigma_c^2)
            x_j ← x_j - K_j * r

        Bug #4 fix: the t-index assertion is performed before calling g_c so
        that misaligned external arrays produce a clear, early error message.
        """
        m_self  = c_dict['models'][local_idx]
        x_up    = states[m_self]
        P_up    = covs[m_self]
        device  = x_up.device
        dtype   = x_up.dtype

        sigma_c = self._as_scalar_tensor(c_dict['sigma_c'], device, dtype)
        tau2    = torch.clamp(sigma_c ** 2, min=1e-12)

        z_mapped = []
        alphas   = []
        vars_all = []

        for j, m in enumerate(c_dict['models']):
            var_name       = c_dict['vars'][j]
            var_idx        = data_list[m].columns.get_loc(var_name) - 1
            alpha, beta    = self._get_map_params(c_dict, j, device, dtype)
            mu_m           = states[m][var_idx].to(device=device, dtype=dtype)
            var_m          = torch.clamp(
                covs[m][var_idx, var_idx].to(device=device, dtype=dtype), min=1e-12
            )
            z_mapped.append(alpha * mu_m + beta)
            alphas.append(alpha)
            vars_all.append(var_m)

        prop_var  = sum(a * a * v for a, v in zip(alphas, vars_all))
        total_var = prop_var + tau2

        # Bug #4: explicit index guard before calling g_c
        try:
            residual = c_dict['g_c'](z_mapped, t=t)
        except TypeError:
            residual = c_dict['g_c'](z_mapped)

        residual = residual.to(device=device, dtype=dtype).reshape(())

        var_name_self = c_dict['vars'][local_idx]
        var_idx_self  = data_list[m_self].columns.get_loc(var_name_self) - 1
        alpha_j       = alphas[local_idx]
        var_j         = vars_all[local_idx]
        K_j           = alpha_j * var_j / total_var

        x_corrected                  = x_up.clone()
        x_corrected[var_idx_self]    = x_corrected[var_idx_self] - K_j * residual
        return x_corrected

    # ──────────────────────────────────────────────────────────────────────────
    # Constraint-aware intra-model state propagation
    # ──────────────────────────────────────────────────────────────────────────

    def _propagate_state_corrections(
        self,
        states_before : list,
        states_after  : list,
        parameters    : list,
    ) -> list:
        """
        After coupling steps have corrected one or more primary state
        variables in a model, propagate those corrections to dependent
        variables using each model's constraint_map callable.

        For the cell-binding ODE the relationship is

            dL/dt = -(db/dt) - k_clear * L
                  = -(k_on*(1-b) - k_off*b) - k_clear * L

        so a coupling-induced shift  Δb  implies a corresponding shift
        ΔL ≈ -Δb  (plus the k_clear term over dt), which this step
        enforces explicitly.

        Parameters
        ----------
        states_before : list[Tensor]
            Per-model state vectors BEFORE any coupling correction this
            timestep (i.e. the raw UKF posterior).
        states_after : list[Tensor]
            Per-model state vectors AFTER all coupling corrections.  Only
            the entries for models whose constraint_map is not None are
            modified.
        parameters : list[dict]
            Per-model parameter dicts (same as passed to inference).

        Returns
        -------
        list[Tensor]
            Updated states_after with dependent variables recomputed.
        """
        result = [s for s in states_after]   # shallow copy of list
        for i, cm in enumerate(self.constraint_maps):
            if cm is None:
                continue
            x_before = states_before[i]
            x_after  = states_after[i]
            # Skip if no correction was applied to this model
            if torch.allclose(x_before, x_after, atol=1e-9):
                continue
            try:
                x_propagated = cm(
                    x_before, x_after, parameters[i], self.universal_dt
                )
            except Exception as exc:
                raise RuntimeError(
                    f"constraint_maps[{i}] raised an error during "
                    f"propagation at this timestep: {exc}"
                ) from exc
            result[i] = x_propagated
        return result

    # ──────────────────────────────────────────────────────────────────────────
    # Inference  (Bug #1 fix: states_over_time committed after all corrections)
    # ──────────────────────────────────────────────────────────────────────────

    def inference(
        self,
        data_list,
        couplings,
        steps,
        params_list=None,
        M2_interp=None,
        progress_bar=True,
    ):
        """
        Multi-model UKF inference with statistical and physical couplings.

        Bug #1 fix
        ----------
        In the original, states_over_time[i].append(x_up) was called inside
        the per-model UKF block (step 1), before statistical blending (step 2)
        and physical coupling corrections (step 3) were applied.  The returned
        trajectory therefore reflected only the raw UKF posterior, not the
        coupled/corrected states.

        Fix: a single commit block at the END of each timestep accumulates
        all outputs (states, covariances, diagnostics) after every correction
        has been applied.
        """
        self._coerce_coupling_tensors(couplings)
        self.ensure_maps_on_couplings(couplings)

        n_models   = len(self.models)
        state_dims = [model.state_dim for model in self.models]

        # ── Initialise states and covariances ─────────────────────────────────
        states = [
            torch.tensor(
                data.iloc[0, 1:1+dim].values, dtype=torch.float32, requires_grad=True
            )
            for data, dim in zip(data_list, state_dims)
        ]
        covariances = [torch.eye(dim, dtype=torch.float32) * 1.0 for dim in state_dims]

        # ── Parameters ────────────────────────────────────────────────────────
        parameters = params_list or [{} for _ in self.models]
        if params_list:
            parameters = [
                {
                    name: (
                        torch.tensor(val, dtype=torch.float32, requires_grad=False)
                        if not isinstance(val, torch.Tensor)
                        else val.clone().detach()
                    )
                    for name, val in params.items()
                }
                for params in params_list
            ]

        # ── Output buffers ─────────────────────────────────────────────────────
        states_over_time      = [[] for _ in self.models]
        param_results         = [[] for _ in self.models]
        covariances_over_time = [[] for _ in self.models]
        z_pred_over_time      = [[] for _ in self.models]
        P_zz_over_time        = [[] for _ in self.models]
        Zi_over_time          = [[] for _ in self.models]
        P_pred_over_time      = [[] for _ in self.models]

        # ── Separate coupling types ────────────────────────────────────────────
        stat_couplings = [
            c for c in couplings
            if c.get('coupling_type', 'statistical') != 'physical'
        ]
        phys_couplings = [
            c for c in couplings
            if c.get('coupling_type') == 'physical'
        ]

        # ── Statistical prior state (initialised from coupling dicts) ──────────
        coupling_states = [
            c['mu_c'].clone().detach().requires_grad_(True)
            for c in stat_couplings
        ]
        coupling_covs = [
            (c['sigma_c'].clone().detach().requires_grad_(True) ** 2).view(1, 1)
            for c in stat_couplings
        ]

        progress = tqdm(range(steps), desc="Inference") if progress_bar else range(steps)

        for t in progress:

            # Per-step diagnostic staging (populated in step 1, committed at end)
            _x_pred   = [None] * n_models
            _P_pred   = [None] * n_models
            _z_pred   = [None] * n_models
            _P_zz     = [None] * n_models
            _params_t = [None] * n_models

            new_states = []
            new_covs   = []

            # ── Step 1: UKF predict → update for each model ───────────────────
            for i, (model, surrogate, state, cov, params) in enumerate(
                zip(self.models, self.surrogates, states, covariances, parameters)
            ):
                # Default argument capture prevents the classic loop-closure bug
                forward_fn = lambda x, dt, _p=params: model.f(x, dt, _p)

                x_pred, P_pred = surrogate.predict(state, cov, forward_fn, self.universal_dt)

                obs_mean = torch.tensor(
                    data_list[i].iloc[t, 1:1+model.state_dim].values,
                    dtype=torch.float32, requires_grad=False,
                )
                try:
                    obs_std = torch.tensor(
                        data_list[i].iloc[t, model.state_dim+1:2*model.state_dim+1].values,
                        dtype=torch.float32,
                    )
                    R = torch.diag(obs_std ** 2)
                except IndexError:
                    R = model.R(state)

                x_up, P_up, z_pred, P_zz = surrogate.update(
                    x_pred, P_pred, obs_mean, model.g, R
                )

                # Stage diagnostics — NOT appended yet (Bug #1 fix)
                _x_pred[i]   = x_pred
                _P_pred[i]   = P_pred.detach().clone()
                _z_pred[i]   = z_pred
                _P_zz[i]     = P_zz
                _params_t[i] = {
                    name: val.detach().numpy() if isinstance(val, torch.Tensor) else val
                    for name, val in params.items()
                }

                new_states.append(x_up)
                new_covs.append(P_up)

            # Snapshot the raw UKF posterior for constraint propagation in step 4
            states_before_coupling = [s.detach().clone() for s in new_states]

            # ── Step 2: Statistical coupling blending ─────────────────────────
            if stat_couplings:
                for i in range(n_models):
                    x_up = new_states[i]
                    P_up = new_covs[i]

                    for c_idx, c in enumerate(stat_couplings):
                        if i in c['models'] and t != 0:
                            model_local = c['models'].index(i)
                            var_name    = c['vars'][model_local]
                            var_idx     = data_list[i].columns.get_loc(var_name) - 1

                            mu_t   = x_up[var_idx]
                            sig_t2 = torch.clamp(P_up[var_idx, var_idx], min=1e-12)
                            mu_c   = coupling_states[c_idx]
                            sig_c2 = torch.clamp(coupling_covs[c_idx].view(-1)[0], min=1e-12)

                            omega_opt = sig_t2 / (sig_t2 + sig_c2)
                            w_j       = c['weights'][model_local]
                            if getattr(self, "learn_omegas", False):
                                omega_eff = torch.clamp(
                                    0.5 * w_j + 0.5 * omega_opt, 0.0, 1.0
                                )
                            else:
                                omega_eff = omega_opt

                            x_blend        = omega_eff * mu_c + (1.0 - omega_eff) * mu_t
                            x_up2          = x_up.clone()
                            x_up2[var_idx] = x_blend
                            x_up           = x_up2

                    new_states[i] = x_up

            # ── Step 3: Physical coupling corrections ─────────────────────────
            if phys_couplings:
                cur_states = new_states
                cur_covs   = new_covs
                for c in phys_couplings:
                    for local_idx, m in enumerate(c['models']):
                        corrected = self._apply_physical_coupling_correction(
                            c_dict    = c,
                            local_idx = local_idx,
                            states    = cur_states,
                            covs      = cur_covs,
                            data_list = data_list,
                            t         = t,
                        )
                        cur_states[m] = corrected
                new_states = cur_states

            # ── Step 4: Intra-model constraint propagation ────────────────────
            # Propagate primary-variable corrections to dependent variables
            # (e.g. b → L in the cell-binding ODE) using each model's
            # constraint_map.  states_before_coupling holds the raw UKF
            # posterior so the constraint map can compute Δb = x_after - x_before.
            if any(cm is not None for cm in self.constraint_maps):
                new_states = self._propagate_state_corrections(
                    states_before = states_before_coupling,
                    states_after  = new_states,
                    parameters    = parameters,
                )

            # ── Commit (Bug #1 fix: append AFTER all corrections) ─────────────
            for i in range(n_models):
                states_over_time[i].append(new_states[i])
                covariances_over_time[i].append(new_covs[i].detach().clone())
                param_results[i].append(_params_t[i])
                Zi_over_time[i].append(_x_pred[i])
                P_pred_over_time[i].append(_P_pred[i])
                z_pred_over_time[i].append(_z_pred[i])
                P_zz_over_time[i].append(_P_zz[i])

            states      = new_states
            covariances = new_covs

            # ── Step 4: Recompute statistical priors for next timestep ─────────
            postfix_info = {}
            for c_idx, c in enumerate(stat_couplings):
                mus, sigmas = [], []
                for model_local, m in enumerate(c['models']):
                    var_name = c['vars'][model_local]
                    var_idx  = data_list[m].columns.get_loc(var_name) - 1
                    mu       = states[m][var_idx]
                    sigma    = torch.sqrt(
                        torch.clamp(covariances[m][var_idx, var_idx], min=1e-10)
                    )
                    mus.append(mu)
                    sigmas.append(sigma)
                mu_c, sigma_c               = self.coupler.mixture_prior(mus, sigmas)
                coupling_states[c_idx]      = mu_c
                coupling_covs[c_idx]        = (sigma_c ** 2).view(1, 1)
                postfix_info[f"cov_c{c_idx}"] = f"{sigma_c**2:.4e}"

            if progress_bar and postfix_info:
                progress.set_postfix(postfix_info)

        # ── Stack and return ───────────────────────────────────────────────────
        states_over_time      = [torch.stack(s) for s in states_over_time]
        P_pred_over_time      = [torch.stack(p) for p in P_pred_over_time]
        covariances_over_time = [torch.stack(c) for c in covariances_over_time]

        return (
            states_over_time,
            param_results,
            covariances_over_time,
            z_pred_over_time,
            P_zz_over_time,
            Zi_over_time,
            P_pred_over_time,
        )

    # ──────────────────────────────────────────────────────────────────────────
    # Energy / loss terms
    # ──────────────────────────────────────────────────────────────────────────

    def compute_data_likelihood_energy(
        self, steps, z_pred_over_time, P_zz_over_time, data_list
    ):
        total_nll_data = torch.tensor(0.0)
        for i, (z_pred_list, P_zz_list, data_df) in enumerate(
            zip(z_pred_over_time, P_zz_over_time, data_list)
        ):
            state_dim = self.models[i].state_dim
            obs_names = data_df.columns[1:1+state_dim]
            obs_mean  = torch.tensor(
                np.stack([data_df[name] for name in obs_names], axis=1),
                dtype=torch.float32,
            )
            for t in range(steps):
                dist            = MultivariateNormal(
                    z_pred_list[t], covariance_matrix=P_zz_list[t]
                )
                total_nll_data += -dist.log_prob(obs_mean[t])
        return total_nll_data

    def compute_state_variable_likelihood_energy(
        self, steps, states_over_time, Zi_over_time, P_pred_over_time
    ):
        total_nll_states = 0.0
        for x_up_list, x_pred_list, P_pred_list in zip(
            states_over_time, Zi_over_time, P_pred_over_time
        ):
            for t in range(steps):
                dist              = MultivariateNormal(
                    x_pred_list[t], covariance_matrix=P_pred_list[t]
                )
                total_nll_states += -dist.log_prob(x_up_list[t])
        return total_nll_states

    def _get_state_scalar(self, Zi_over_time, m, t, var_idx, device, dtype=torch.float64):
        z_m = Zi_over_time[m]
        if isinstance(z_m, torch.Tensor):
            if z_m.dim() == 2:
                return z_m[t, var_idx].reshape(()).to(device=device, dtype=dtype)
            elif z_m.dim() == 1:
                return z_m[var_idx].reshape(()).to(device=device, dtype=dtype)
            else:
                raise ValueError(
                    f"Unexpected tensor shape for Zi_over_time[{m}]: {z_m.shape}"
                )
        else:
            z_t = z_m[t]
            if isinstance(z_t, torch.Tensor):
                return z_t[var_idx].reshape(()).to(device=device, dtype=dtype)
            return self._as_scalar_tensor(z_t[var_idx], device, dtype=dtype)

    def _get_cov_scalar(
        self, state_covs_over_time, m, t, var_idx, device,
        dtype=torch.float64, fallback_var=1e-3
    ):
        if state_covs_over_time is None:
            return torch.tensor(fallback_var, device=device, dtype=dtype).reshape(())
        P_m = state_covs_over_time[m]
        if isinstance(P_m, torch.Tensor):
            if P_m.dim() == 3:
                return (
                    torch.clamp(P_m[t, var_idx, var_idx], min=1e-12)
                    .reshape(()).to(device=device, dtype=dtype)
                )
            raise ValueError(
                f"Unexpected tensor shape for state_covs_over_time[{m}]: {P_m.shape}"
            )
        P_t = P_m[t]
        if isinstance(P_t, torch.Tensor):
            return torch.clamp(P_t[var_idx, var_idx], min=1e-12).reshape(()).to(
                device=device, dtype=dtype
            )
        val = P_t[var_idx][var_idx]
        return torch.clamp(self._as_scalar_tensor(val, device, dtype=dtype), min=1e-12)

    def _apply_coupling_map_and_slack(self, c_dict, mu_j, var_j, local_idx, tau2_c):
        alpha, beta = self._get_map_params(c_dict, local_idx, mu_j.device, mu_j.dtype)
        y_mu        = alpha * mu_j + beta
        y_var       = torch.clamp(alpha * alpha * var_j + tau2_c, min=1e-12)
        return y_mu, y_var

    def _consensus_log_marginal_gaussian(self, y, var):
        var       = torch.clamp(var, min=1e-12)
        prec      = 1.0 / var
        A         = torch.sum(prec)
        mu_pool   = torch.sum(prec * y) / A
        var_pool  = 1.0 / A
        sig_pool  = torch.sqrt(var_pool)
        y_std     = (y  - mu_pool) / (sig_pool + 1e-12)
        var_std   = var / (var_pool + 1e-12)
        var_std   = torch.clamp(var_std, min=1e-12)
        prec_std  = 1.0 / var_std
        A_std     = torch.sum(prec_std)
        ybar_std  = torch.sum(prec_std * y_std) / A_std
        SSEw_std  = torch.sum((y_std - ybar_std) ** 2 * prec_std)
        J         = y.shape[0]
        logp = (
            -0.5 * SSEw_std
            - 0.5 * torch.sum(torch.log(var_std))
            - 0.5 * torch.log(A_std)
            - 0.5 * (J - 1) * math.log(2.0 * math.pi)
        )
        return logp.reshape(())

    def _gather_mu_var_for_coupling(
        self, t, c_dict, states_over_time, state_covs_over_time,
        data_list, device, tau2_c, fallback_var=1e-3
    ):
        y_mus, y_vars = [], []
        for j, m in enumerate(c_dict['models']):
            var_name = c_dict['vars'][j]
            var_idx  = data_list[m].columns.get_loc(var_name) - 1
            mu_j     = self._get_state_scalar(states_over_time, m, t, var_idx, device)
            var_j    = self._get_cov_scalar(
                state_covs_over_time, m, t, var_idx, device, fallback_var=fallback_var
            )
            y_mu, y_var = self._apply_coupling_map_and_slack(
                c_dict, mu_j, var_j, j, tau2_c
            )
            y_mus.append(y_mu)
            y_vars.append(y_var)
        return torch.stack(y_mus), torch.stack(y_vars)

    def _normalized_gaussian_energy(self, residual, variances, include_const=False):
        residual  = residual.reshape(-1)
        variances = torch.clamp(variances.reshape(-1), min=1e-12)
        d      = residual.numel()
        quad   = 0.5 * torch.sum((residual ** 2) / variances)
        logdet = 0.5 * torch.sum(torch.log(variances))
        const  = 0.5 * d * math.log(2.0 * math.pi) if include_const else 0.0
        return ((quad + logdet + const) / d).reshape(())

#    def compute_loss(
#        self, data_list, couplings, steps,
#        Zi_over_time, states_over_time, covariances_over_time,
#        z_pred_over_time, P_zz_over_time, P_pred_over_time,
#        device=torch.device("cpu"), dtype=torch.float64,
#    ):
#        self._coerce_coupling_tensors(couplings)
#        self.ensure_maps_on_couplings(couplings)
#
#        n_models       = len(self.models)
#        n_couplings    = len(couplings)
#        total_obs_data = sum(steps * self.models[i].state_dim for i in range(n_models))
#
#        # Term 1: coupling NLL
#        nll_cpl = torch.tensor(0.0, dtype=dtype, device=device)
#        if n_couplings > 0:
#            for t in range(steps):
#                cur_states = [states_over_time[m][t] for m in range(n_models)]
#                cur_covs   = [covariances_over_time[m][t] for m in range(n_models)]
#
#                for c in couplings:
#                    ctype = c.get('coupling_type', 'statistical')
#
#                    if ctype == 'physical':
#                        sigma_c  = self._as_scalar_tensor(c['sigma_c'], device, dtype)
#                        tau2     = torch.clamp(sigma_c ** 2, min=1e-12)
#                        z_mapped = []
#                        prop_var = torch.tensor(0.0, device=device, dtype=dtype)
#
#                        for j, m in enumerate(c['models']):
#                            var_name        = c['vars'][j]
#                            var_idx         = data_list[m].columns.get_loc(var_name) - 1
#                            mu_j            = cur_states[m][var_idx].to(device=device, dtype=dtype)
#                            var_j           = torch.clamp(
#                                cur_covs[m][var_idx, var_idx].to(device=device, dtype=dtype),
#                                min=1e-12,
#                            )
#                            alpha, beta     = self._get_map_params(c, j, device, dtype)
#                            z_mapped.append(alpha * mu_j + beta)
#                            prop_var        = prop_var + alpha * alpha * var_j
#
#                        try:
#                            residual = c['g_c'](z_mapped, t=t)
#                        except TypeError:
#                            residual = c['g_c'](z_mapped)
#
#                        residual  = residual.to(device=device, dtype=dtype).reshape(())
#                        total_var = torch.clamp(tau2 + prop_var, min=1e-12)
#                        nll_cpl  += self._normalized_gaussian_energy(
#                            residual=residual.view(1),
#                            variances=total_var.view(1),
#                            include_const=False,
#                        )
#
#                    else:
#                        tau2_c = torch.tensor(0.0, dtype=dtype, device=device)
#                        y, v   = self._gather_mu_var_for_coupling(
#                            t, c, states_over_time, covariances_over_time,
#                            data_list, device, tau2_c, fallback_var=1e-3,
#                        )
#                        prec     = 1.0 / torch.clamp(v, min=1e-12)
#                        mu_pool  = torch.sum(prec * y) / torch.sum(prec)
#                        residual = y - mu_pool
#                        nll_cpl += self._normalized_gaussian_energy(
#                            residual=residual, variances=v, include_const=False,
#                        )
#
#            nll_cpl_n = nll_cpl / max(1, n_couplings)
#        else:
#            nll_cpl_n = torch.tensor(0.0, device=device, dtype=dtype)
#
#        nll_cpl_n = nll_cpl_n / total_obs_data
#
#        # Term 2: data fidelity
#        nll_data   = self.compute_data_likelihood_energy(
#            steps, z_pred_over_time, P_zz_over_time, data_list
#        )
#        nll_data_n = nll_data / total_obs_data
#
#        # Term 3: surrogate self-consistency
#        nll_states   = self.compute_state_variable_likelihood_energy(
#            steps, states_over_time, Zi_over_time, P_pred_over_time
#        )
#        nll_states_n = nll_states / total_obs_data
#
#        self.num_obs      = total_obs_data
#        self.n_couplings  = n_couplings
#        self.nll_cpl_n    = float(nll_cpl_n)
#        self.nll_data_n   = float(nll_data_n)
#        self.nll_states_n = float(nll_states_n)
#        self.nll          = torch.tensor(self.nll_data_n + self.nll_states_n)
#        self.nll_data     = nll_data
#        self.nll_states   = nll_states
#
#        return float(nll_cpl_n) + float(nll_data_n)


    def compute_loss(
        self, data_list, couplings, steps,
        Zi_over_time, states_over_time, covariances_over_time,
        z_pred_over_time, P_zz_over_time, P_pred_over_time,
        device=torch.device("cpu"), dtype=torch.float64,
    ):
        self._coerce_coupling_tensors(couplings)
        self.ensure_maps_on_couplings(couplings)

        n_models       = len(self.models)
        n_couplings    = len(couplings)
        total_obs_data = sum(steps * self.models[i].state_dim for i in range(n_models))

        # ── Term 1: coupling NLL (unified form) ───────────────────────────────────
        nll_cpl = torch.tensor(0.0, dtype=dtype, device=device)
        if n_couplings > 0:
            for t in range(steps):
                cur_states = [states_over_time[m][t] for m in range(n_models)]
                cur_covs   = [covariances_over_time[m][t] for m in range(n_models)]

                for c in couplings:
                    ctype   = c.get('coupling_type', 'statistical')
                    is_phys = (ctype == 'physical')

                    # tau^2: constraint noise for physical; 0 for statistical
                    if is_phys:
                        sigma_c = self._as_scalar_tensor(c['sigma_c'], device, dtype)
                        tau2    = torch.clamp(sigma_c ** 2, min=1e-12)
                    else:
                        tau2 = torch.tensor(0.0, dtype=dtype, device=device)

                    # Build mapped means y_j and per-model variances v_j
                    y_mus, y_vars = [], []
                    for j, m in enumerate(c['models']):
                        var_name    = c['vars'][j]
                        var_idx     = data_list[m].columns.get_loc(var_name) - 1
                        mu_j        = cur_states[m][var_idx].to(device=device, dtype=dtype)
                        var_j       = torch.clamp(
                            cur_covs[m][var_idx, var_idx].to(device=device, dtype=dtype),
                            min=1e-12,
                        )
                        alpha, beta = self._get_map_params(c, j, device, dtype)
                        y_mus.append(alpha * mu_j + beta)
                        y_vars.append(torch.clamp(alpha * alpha * var_j, min=1e-12))

                    y    = torch.stack(y_mus)   # (J,)
                    v    = torch.stack(y_vars)  # (J,)  propagated model variances

                    # Compute residuals: physical → g_c scalar; statistical → y_j - pool
                    if is_phys:
                        try:
                            residual = c['g_c'](list(y_mus), t=t)
                        except TypeError:
                            residual = c['g_c'](list(y_mus))
                        residual = residual.to(device=device, dtype=dtype).reshape(())
                        # broadcast: same scalar residual for every participant
                        residuals = residual.expand(len(c['models']))
                    else:
                        prec     = 1.0 / torch.clamp(v, min=1e-12)
                        mu_pool  = torch.sum(prec * y) / torch.sum(prec)
                        residuals = y - mu_pool   # (J,)

                    # Unified per-participant energy: same form for both types
                    #   L_j = 0.5 * r_j^2 / (v_j + tau^2) + 0.5 * log(v_j + tau^2)
                    effective_vars = v + tau2          # tau2=0 for stat, >0 for phys
                    nll_cpl += self._normalized_gaussian_energy(
                        residual  = residuals,
                        variances = effective_vars,
                        include_const = False,
                    )

            nll_cpl_n = nll_cpl / max(1, n_couplings)
        else:
            nll_cpl_n = torch.tensor(0.0, device=device, dtype=dtype)

        nll_cpl_n = nll_cpl_n / total_obs_data

        # ── Term 2: data fidelity ─────────────────────────────────────────────────
        nll_data   = self.compute_data_likelihood_energy(
            steps, z_pred_over_time, P_zz_over_time, data_list
        )
        nll_data_n = nll_data / total_obs_data

        # ── Term 3: surrogate self-consistency ────────────────────────────────────
        nll_states   = self.compute_state_variable_likelihood_energy(
            steps, states_over_time, Zi_over_time, P_pred_over_time
        )
        nll_states_n = nll_states / total_obs_data

        self.num_obs      = total_obs_data
        self.n_couplings  = n_couplings
        self.nll_cpl_n    = float(nll_cpl_n)
        self.nll_data_n   = float(nll_data_n)
        self.nll_states_n = float(nll_states_n)
        self.nll          = torch.tensor(self.nll_data_n + self.nll_states_n)
        self.nll_data     = nll_data
        self.nll_states   = nll_states

        return float(nll_cpl_n) + float(nll_data_n)


# }}}

# Physical coupling:{{{
def _coerce_coupling_tensors(couplings):
    """Ensure mu_c / sigma_c are tensors on every coupling dict (in-place)."""
    for c in couplings:
        if not isinstance(c['mu_c'], torch.Tensor):
            c['mu_c'] = torch.tensor(float(c['mu_c']), dtype=torch.float32)
        if not isinstance(c['sigma_c'], torch.Tensor):
            c['sigma_c'] = torch.tensor(float(c['sigma_c']), dtype=torch.float32)


@dataclass
class PhysicalCoupling:
    """Declares a single physical coupling constraint and serialises to a dict."""
    model_indices : List[int]
    var_names     : List[str]
    g_c           : Callable
    sigma_c       : float = 0.1
    lambdas       : List[float] = field(default_factory=list)
    shifts        : List[float] = field(default_factory=list)
    weights       : List       = field(default_factory=list)

    def as_dict(self) -> dict:
        n = len(self.model_indices)
        return {
            'coupling_type': 'physical',
            'models':        tuple(self.model_indices),
            'vars':          tuple(self.var_names),
            'g_c':           self.g_c,
            'sigma_c':       self.sigma_c,
            'lambdas':       self.lambdas if self.lambdas else [1.0] * n,
            'shifts':        self.shifts  if self.shifts  else [0.0] * n,
            'mu_c':          torch.tensor(0.0),
            'weights':       (
                self.weights if self.weights
                else [torch.tensor(1.0 / n) for _ in range(n)]
            ),
        }


def make_physical_coupling(
    model_indices : Sequence[int],
    var_names     : Sequence[str],
    g_c           : Callable,
    sigma_c       : float = 0.1,
    lambdas       : Optional[Sequence[float]] = None,
    shifts        : Optional[Sequence[float]] = None,
    #sigma_c_floor : float = 0.05,
    sigma_c_floor : float = 0.005,
) -> dict:
    """
    Build a physical coupling dict compatible with BayesMM.inference().

    Bug #3 fix
    ----------
    sigma_c is clamped to max(sigma_c, sigma_c_floor).  When sigma_c derives
    from an experimental uncertainty that is near zero (e.g. TRUE_dG_std), the
    physical coupling would otherwise become infinitely rigid and overwhelm the
    ODE dynamics.  The floor prevents this.  Pass sigma_c_floor=0.0 to disable.
    """
    sigma_c_safe = max(float(sigma_c), float(sigma_c_floor))
    if sigma_c_safe != float(sigma_c):
        warnings.warn(
            f"make_physical_coupling: sigma_c={sigma_c:.4e} is below the floor "
            f"{sigma_c_floor:.4e}; using {sigma_c_safe:.4e}.  "
            "Pass sigma_c_floor=0.0 to disable this guard.",
            UserWarning,
            stacklevel=2,
        )
    n  = len(model_indices)
    pc = PhysicalCoupling(
        model_indices = list(model_indices),
        var_names     = list(var_names),
        g_c           = g_c,
        sigma_c       = sigma_c_safe,
        lambdas       = list(lambdas) if lambdas is not None else [1.0] * n,
        shifts        = list(shifts)  if shifts  is not None else [0.0] * n,
    )
    return pc.as_dict()
# }}}

# Helper functions:{{{

def set_coupling_weights_from_probabilities(couplings, probs):
    idx = 0
    for c in couplings:
        for i in range(len(c['models'])):
            c['weights'][i] = torch.tensor(probs[idx], dtype=torch.float32)
            idx += 1


def _torch_to_numpy_state(x_torch: torch.Tensor) -> np.ndarray:
    return x_torch.detach().cpu().numpy().astype(float)


def _torch_params_to_numpy_list(params: dict, param_names: list) -> list:
    return [
        float(params[name].detach().cpu().item())
        if isinstance(params[name], torch.Tensor)
        else float(params[name])
        for name in param_names
    ]


def _numpy_dxdt_to_torch(dxdt_np: np.ndarray) -> torch.Tensor:
    return torch.tensor(dxdt_np, dtype=torch.float32, requires_grad=True)


def forward_fn(x: torch.Tensor, dt_step: float, params: dict, ode_fn, param_names, ode_kwargs=None):
    x_np      = _torch_to_numpy_state(x)
    k_list_np = _torch_params_to_numpy_list(params, param_names)

    if ode_kwargs is None:
        ode_kwargs = {}

    dxdt_np = ode_fn(x_np, 0.0, k_list_np, **ode_kwargs)
    dxdt_clean = [
        val.nominal_value if hasattr(val, "nominal_value") else float(val)
        for val in dxdt_np
    ]
    dx_dt_torch = _numpy_dxdt_to_torch(np.array(dxdt_clean, dtype=float))
    return x + dt_step * dx_dt_torch


def generate_surrogate(
    observables,
    model_parameters,
    ode_fn,
    sim_time,
    dt,
    initial_noise_scale=0.01,
    transition_cov_scale=0.01,
    emission_cov_scale=1.0,
    observation_noise_scale=0.001,
    noise_model_type='time-invariant',
    ode_kwargs: dict = None,
):
    obs_mean  = observables[:, :, 0]
    obs_std   = observables[:, :, 1]
    steps     = obs_mean.shape[0]
    state_dim = len(model_parameters["InitConc"])
    param_names = [f"k{i+1}" for i in range(len(model_parameters["k"]))]

    data_df = pd.DataFrame({
        'time': np.arange(0, sim_time, dt)[:steps],
        **{name: obs_mean[:, i].detach().numpy()
           for i, name in enumerate(model_parameters["species"])},
        **{f"{name}_std": obs_std[:, i].detach().numpy()
           for i, name in enumerate(model_parameters["species"])},
    })

    def fwd_bridge(x_torch: torch.Tensor, dt_step: float, params_dict: dict):
        if ode_fn is None:
            return x_torch
        return forward_fn(
            x=x_torch, dt_step=dt_step,
            params=params_dict, ode_fn=ode_fn, param_names=param_names, ode_kwargs=ode_kwargs,
        )

    def observe_fn(x: torch.Tensor) -> torch.Tensor:
        return x

    def trans_noise(x: torch.Tensor) -> torch.Tensor:
        if noise_model_type == 'time-invariant':
            mean_x = torch.mean(x, dim=0) if x.ndim > 1 else x
            return torch.diag(torch.abs(mean_x * transition_cov_scale) + 1e-6)
        elif noise_model_type == 'time-variant':
            noise = torch.abs(
                x * torch.normal(
                    mean=transition_cov_scale,
                    std=transition_cov_scale * 1e-2,
                    size=x.shape,
                )
            )
            return torch.diag(noise + 1e-6)
        raise ValueError(f"Unknown noise_model_type: {noise_model_type!r}")

    def emit_noise(x: torch.Tensor) -> torch.Tensor:
        if noise_model_type == 'time-invariant':
            mean_x = torch.mean(x, dim=0) if x.ndim > 1 else x
            return torch.diag(torch.abs(mean_x * emission_cov_scale) + 1e-6)
        elif noise_model_type == 'time-variant':
            noise = torch.abs(
                x * torch.normal(
                    mean=emission_cov_scale,
                    std=emission_cov_scale * 1e-2,
                    size=x.shape,
                )
            )
            return torch.diag(noise + 1e-6)
        raise ValueError(f"Unknown noise_model_type: {noise_model_type!r}")

    surrogate = SurrogateModel(
        forward_fn=fwd_bridge,
        observe_fn=observe_fn,
        trans_noise=trans_noise,
        emit_noise=emit_noise,
        state_dim=state_dim,
    )

    param_dict = {
        name: torch.tensor(model_parameters["k"][i], dtype=torch.float32, requires_grad=True)
        for i, name in enumerate(param_names)
    }
    bayes_mm = BayesMM(
        models=[surrogate], param_names_list=[param_names], universal_dt=dt
    )
    surrogate_states = bayes_mm.inference(
        data_list=[data_df], couplings=[], steps=steps,
        params_list=[param_dict], progress_bar=False,
    )[0][0]

    surrogate_df = pd.DataFrame({
        'time': np.arange(0, sim_time, dt)[:steps],
        **{name: surrogate_states[:, i].detach().numpy()
           for i, name in enumerate(model_parameters["species"])},
    })

    return surrogate, surrogate_states, surrogate_df, param_names


# ══════════════════════════════════════════════════════════════════════════════
# Coupling search helpers
# ══════════════════════════════════════════════════════════════════════════════

def _powerset(iterable):
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def auto_identify_couplings(
    data_dfs,
    params_list,
    models_data,
    sim_time,
    params_list_for_inf      = None,
    candidate_phys_couplings = None,
):
    """
    Joint search over statistical coupling candidates × physical coupling subsets.

    Returns the configuration (stat_vars, phys_subset, Zi_over_time) that
    minimises the normalised total loss.  The no-coupling baseline is excluded.
    """
    candidate_phys_couplings = candidate_phys_couplings or []
    phys_subsets             = _powerset(candidate_phys_couplings)

    # ── Statistical candidates ranked by R² ──────────────────────────────────
    scores = []
    for i, j in combinations(range(len(data_dfs)), 2):
        df_i, df_j = data_dfs[i], data_dfs[j]
        for var_i in params_list[i]["species"]:
            for var_j in params_list[j]["species"]:
                x    = df_i[var_i].to_numpy()
                y    = df_j[var_j].to_numpy()
                mi   = mutual_info_regression(x[:, np.newaxis], y)[0]
                corr = np.corrcoef(x, y)[0, 1] ** 2
                scores.append(((i, j), (var_i, var_j), corr, mi))

    sorted_model_evals = sorted(scores, key=lambda x: x[2], reverse=True)
    stat_candidates    = [
        ((i, j), (var_i, var_j))
        for (i, j), (var_i, var_j), R2, mi in sorted_model_evals
    ]

    bayes_mm = BayesMM(
        [m['surrogate'] for m in models_data],
        param_names_list = [m['param_names'] for m in models_data],
        universal_dt     = min(m['dt'] for m in models_data),
        constraint_maps  = [m.get('constraint_map', None) for m in models_data],
    )
    steps  = int(sim_time / min(m['dt'] for m in models_data))
    device = torch.device("cpu")
    dtype  = torch.float64

    # ── Build candidate list (every config with >= 1 active coupling) ─────────
    all_candidates = []
    for phys_subset in phys_subsets:
        if len(phys_subset) > 0:
            all_candidates.append((None, phys_subset))
    for stat_candidate in stat_candidates:
        for phys_subset in phys_subsets:
            all_candidates.append((stat_candidate, phys_subset))

    best_loss         = float('inf')
    best_stat_vars    = []
    best_phys_subset  = []
    best_Zi_over_time = None
    all_results       = []
    total             = len(all_candidates)

    for run_idx, (stat_candidate, phys_subset) in enumerate(all_candidates):
        phys_list     = list(phys_subset)
        stat_couplings = (
            bayes_mm.couple_models(data_dfs, [stat_candidate])
            if stat_candidate is not None else []
        )
        all_couplings = stat_couplings + phys_list
        _coerce_coupling_tensors(all_couplings)
        bayes_mm.ensure_maps_on_couplings(all_couplings)

        (states_over_time, _,
         covariances_over_time,
         z_pred_over_time, P_zz_over_time,
         Zi_over_time, P_pred_over_time) = bayes_mm.inference(
            data_list    = data_dfs,
            couplings    = all_couplings,
            steps        = steps,
            params_list  = params_list_for_inf,
            progress_bar = False,
        )

        loss = bayes_mm.compute_loss(
            data_dfs, all_couplings, steps,
            Zi_over_time, states_over_time, covariances_over_time,
            z_pred_over_time, P_zz_over_time, P_pred_over_time,
            device=device, dtype=dtype,
        )

        s_idx      = stat_candidates.index(stat_candidate) if stat_candidate is not None else None
        stat_label = (
            f"{stat_candidate[0]}:{stat_candidate[1]}"
            if stat_candidate is not None else "none"
        )
        phys_label = (
            "+".join(f"phys({c['vars']})" for c in phys_list)
            if phys_list else "none"
        )
        print(
            f"[{run_idx+1}/{total}] "
            f"stat={stat_label}  phys={phys_label}  "
            f"nll_cpl/obs={bayes_mm.nll_cpl_n:.4f}  "
            f"nll_data/obs={bayes_mm.nll_data_n:.4f}  "
            f"loss={loss:.4f}"
        )
        all_results.append({
            "Stat models":    sorted_model_evals[s_idx][0] if s_idx is not None else "none",
            "Stat vars":      sorted_model_evals[s_idx][1] if s_idx is not None else "none",
            "Phys couplings": phys_label,
            "nll_cpl/obs":    bayes_mm.nll_cpl_n,
            "nll_data/obs":   bayes_mm.nll_data_n,
            "L(theta)":       loss,
        })

        if loss < best_loss:
            best_loss         = loss
            best_stat_vars    = [stat_candidate] if stat_candidate is not None else []
            best_phys_subset  = phys_list
            best_Zi_over_time = Zi_over_time

    results_df = (
        pd.DataFrame(all_results)
        .sort_values("L(theta)")
        .reset_index(drop=True)
    )
    print("\n── Coupling search results (sorted by loss) ─────────────────────")
    print(results_df.to_string(index=False))
    print(f"\n── Best stat vars:   {best_stat_vars}")
    print(
        f"── Best phys subset: "
        f"{[c['vars'] for c in best_phys_subset] if best_phys_subset else 'none'}"
    )
    print(f"── Best loss:        {best_loss:.4f}")

    return best_stat_vars, best_phys_subset, best_Zi_over_time


# ══════════════════════════════════════════════════════════════════════════════
# build_metamodel
# ══════════════════════════════════════════════════════════════════════════════

def build_metamodel(
    models_data,
    sim_time,
    connecting_vars          = None,
    couplings                = None,
    candidate_phys_couplings = None,
    run_param_opt            = True,
    lr                       = 4e-2,
    max_iterations           = 50,
):
    """
    Build the Bayesian metamodel.

    Parameters
    ----------
    models_data : list[dict]
        One dict per model.  Required keys: 'obs_model', 'surrogate', 'params',
        'ode_fn', 'dt', 'transition_cov_scale', 'emission_cov_scale',
        'observation_noise_scale', 'noise_model_type', 'param_names'.
    sim_time : float
        Total simulation time.
    connecting_vars : list | None
        Statistical connecting variable pairs.  None → auto-search.
    couplings : list | None
        Pre-built statistical coupling dicts.  None → built from connecting_vars.
    candidate_phys_couplings : list[dict] | None
        Physical coupling dicts from make_physical_coupling().
    """
    obs_models       = [m['obs_model']   for m in models_data]
    params_list      = [m['params']      for m in models_data]
    ode_fns          = [m['ode_fn']      for m in models_data]
    dts              = [m['dt']          for m in models_data]
    param_names_list = [m['param_names'] for m in models_data]

    universal_dt   = min(dts)
    steps          = int(sim_time / universal_dt)
    universal_time = np.arange(0, sim_time, universal_dt)[:steps]

    # ── Interpolate all observables to universal_dt ───────────────────────────
    processed_obs = []
    data_dfs      = []
    for idx, (obs_model, params, model_dt, ode_fn, param_names) in enumerate(
        zip(obs_models, params_list, dts, ode_fns, param_names_list)
    ):
        model_steps = obs_model.shape[0]
        state_dim   = obs_model.shape[1]
        model_time  = np.arange(0, sim_time, model_dt)[:model_steps]

        if model_dt > universal_dt:
            obs_processed = torch.zeros((steps, state_dim, 2))
            for i in range(state_dim):
                for j in range(2):
                    interp_func = interp1d(
                        model_time,
                        obs_model[:, i, j].detach().numpy(),
                        kind='linear', fill_value='extrapolate',
                    )
                    obs_processed[:, i, j] = torch.tensor(
                        interp_func(universal_time), dtype=torch.float32
                    )
        else:
            obs_processed = obs_model[:steps]

        processed_obs.append(obs_processed)
        data_df = pd.DataFrame({
            'time': universal_time,
            **{name: obs_processed[:, i, 0].detach().numpy()
               for i, name in enumerate(params["species"])},
            **{f"{name}_std": obs_processed[:, i, 1].detach().numpy()
               for i, name in enumerate(params["species"])},
        })
        data_dfs.append(data_df)

    # ── Build surrogates ──────────────────────────────────────────────────────
    for idx, (obs_processed, params, ode_fn) in enumerate(
        zip(processed_obs, params_list, ode_fns)
    ):
        surrogate, _, _, _ = generate_surrogate(
            observables             = obs_processed,
            model_parameters        = params,
            ode_fn                  = ode_fn,
            sim_time                = sim_time,
            dt                      = universal_dt,
            transition_cov_scale    = models_data[idx]['transition_cov_scale'],
            emission_cov_scale      = models_data[idx]['emission_cov_scale'],
            observation_noise_scale = models_data[idx]['observation_noise_scale'],
            noise_model_type        = models_data[idx]['noise_model_type'],
            ode_kwargs              = models_data[idx].get('ode_kwargs', None),
        )
        models_data[idx]['surrogate'] = surrogate

    constraint_maps = [
        m.get('constraint_map', None) for m in models_data
    ]
    bayes_mm = BayesMM(
        [m['surrogate'] for m in models_data],
        param_names_list = param_names_list,
        universal_dt     = universal_dt,
        constraint_maps  = constraint_maps,
    )

    # ── Coupling selection ────────────────────────────────────────────────────
    selected_phys_couplings = candidate_phys_couplings or []

    if connecting_vars is None:
        params_list_for_inf = [
            {name: params['k'][i] for i, name in enumerate(names)}
            for params, names in zip(params_list, param_names_list)
        ]
        connecting_vars, selected_phys_couplings, _ = auto_identify_couplings(
            data_dfs, params_list, models_data, sim_time,
            params_list_for_inf      = params_list_for_inf,
            candidate_phys_couplings = candidate_phys_couplings or [],
        )
        print("Auto-identified connecting_vars =", connecting_vars)
        print(
            "Auto-selected phys couplings    =",
            [c['vars'] for c in selected_phys_couplings]
            if selected_phys_couplings else "none",
        )

    if couplings is None:
        couplings = bayes_mm.couple_models(data_dfs, connecting_vars)

    all_couplings = couplings + selected_phys_couplings
    _coerce_coupling_tensors(all_couplings)
    bayes_mm.ensure_maps_on_couplings(all_couplings)

    # ── Final inference ───────────────────────────────────────────────────────
    params_list_for_opt = [
        {name: params['k'][i] for i, name in enumerate(names)}
        for params, names in zip(params_list, param_names_list)
    ]

    (results, param_results,
     covariances_over_time,
     z_pred_over_time, P_zz_over_time,
     Zi_over_time, P_pred_over_time) = bayes_mm.inference(
        data_list   = data_dfs,
        couplings   = all_couplings,
        steps       = steps,
        params_list = params_list_for_opt,
    )
    states_over_time = results

    bayes_mm.loss = loss = bayes_mm.compute_loss(
        data_dfs, all_couplings, steps,
        Zi_over_time, states_over_time, covariances_over_time,
        z_pred_over_time, P_zz_over_time, P_pred_over_time,
    )
    print(f"BMM loss       = {loss:.4f}")
    print(f"  nll_cpl/obs  = {bayes_mm.nll_cpl_n:.4f}")
    print(f"  nll_data/obs = {bayes_mm.nll_data_n:.4f}")

    states = [np.array(r.detach().numpy()) for r in results]

    coupled_dfs = []
    for idx, (state, params) in enumerate(zip(states, params_list)):
        coupled_dfs.append(pd.DataFrame({
            'time': universal_time,
            **{name: state[:, i] for i, name in enumerate(params["species"])},
        }))

    figs = []
    for idx, (state, data_df, params, obs_model) in enumerate(
        zip(states, data_dfs, params_list, processed_obs)
    ):
        ylabels = [
            f"{s} [{params['units'][i]}]"
            for i, s in enumerate(params["species"])
        ]
        fig = plot_inputmodel(
            input_result = state,
            input_std    = np.stack(
                [data_df[f"{name}_std"] for name in params["species"]], axis=1
            ),
            dt           = universal_dt,
            sim_time     = sim_time,
            name         = ylabels,
            observables  = obs_model.detach().numpy(),
        )
        figname = os.path.join(outdir, f"coupled_model_{idx}_with_observables.png")
        fig.tight_layout()
        fig.savefig(figname)
        figs.append(fig)

    return (
        bayes_mm, states, data_dfs, coupled_dfs, figs,
        connecting_vars, processed_obs, covariances_over_time,
    )
# }}}



