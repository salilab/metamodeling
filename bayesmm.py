# Import libraries:{{{
from __future__ import annotations
import os, math, re, warnings, time, copy, inspect
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributions as dist
import matplotlib.pyplot as plt
from matplotlib import patches
from tqdm import trange,tqdm
from typing import Optional, Union, Callable, Tuple, List, Dict, Iterable, Any, Sequence
from torchdiffeq import odeint
from torch.optim.lr_scheduler import ReduceLROnPlateau
from matplotlib.pyplot import MultipleLocator
from itertools import combinations
from sklearn.feature_selection import mutual_info_regression
from torch.distributions import MultivariateNormal
import torch.nn.functional as F
from scipy.interpolate import interp1d
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

# }}}

# ODEs: {{{

def get_observations(rate_eqs, parameters, dt, sim_time, noise_scale=0.0, kwargs=None):
    """
    Runs the uncertainty-aware ODE model and returns observations with mean and std.
    Args:
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
    steps = int(np.ceil(sim_time / dt))
    df, _, sol = simulate_model(rate_eqs=rate_eqs, parameters=parameters,
                                steps=steps, simulation_time=sim_time, kwargs=kwargs)
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

    if noise_scale is not None and noise_scale != 0.0:
        obs_noise = np.abs(noise_scale * means_np)
        total_std = np.sqrt(stds_np**2 + obs_noise**2) + 1e-6
    else:
        total_std = stds_np + 1e-6  # ensure strictly positive

    # pack into (steps, n_species, 2)
    obs_mean_tensor = torch.tensor(means_np, dtype=torch.float32)
    obs_std_tensor  = torch.tensor(total_std, dtype=torch.float32)
    observations = torch.stack([obs_mean_tensor, obs_std_tensor], dim=2)
    t_grid = df["time"].to_numpy()

    return observations, t_grid, df


def _to_list(x):
    if isinstance(x, np.ndarray): return [xi for xi in x.tolist()]
    else: return [xi for xi in x]

def _sanitize_state(y_arr):
    cleaned = []
    for val in y_arr:
        # ensure val is ufloat
        if hasattr(val, "nominal_value"):
            nom,std = val.nominal_value, val.std_dev
        else:
            nom,std = float(val), 0.0
            val = ufloat(nom, std)
        #if (not np.isfinite(nom)) or (nom < 0.0):
        #    cleaned.append(ufloat(0.0, 0.0))
        #else:
        #    cleaned.append(val)
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

        if dt <= 0: raise ValueError(f"Non-positive dt encountered at step {idx}: dt={dt}")

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

        y_next = y_curr + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # RK4 update
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

def simulate_model( rate_eqs, parameters, steps, simulation_time, kwargs=None):
    t_grid = np.linspace(0.0, simulation_time, steps)
    # Build initial conditions with uncertainties
    init_conc = []
    for i in range(len(parameters["InitConc"])):
        mean_i = parameters["InitConc"][i]
        std_i  = parameters["InitConc_std"][i]
        if std_i is not None:
            init_conc.append(u.ufloat(mean_i, std_i))
        else:
            init_conc.append(u.ufloat(mean_i, 0.0))
    # rate constants with uncertainties
    k_list = []
    for i in range(len(parameters["params"])):
        mean_k = parameters["params"][i]
        std_k  = parameters["params_std"][i]
        if std_k is not None:
            k_list.append(u.ufloat(mean_k, std_k))
        else:
            k_list.append(u.ufloat(mean_k, 0.0))
    # Resample observable if provided
    if kwargs is None:
        kwargs = {}
    if "obs_fn" in kwargs:
        raw = kwargs["obs_fn"]
        vals = raw(t_grid)
        kwargs["obs_fn"] = interp1d(t_grid, vals, kind="linear", bounds_error=False, fill_value=(vals[0], vals[-1]))

    # Integrate with custom RK4 that keeps ufloats alive
    # NOTE: rk4_integrate must internally clamp negatives and replace NaN/inf,
    # because np.clip/np.nan_to_num won't like dtype=object arrays of ufloats.
    sol = rk4_integrate(rate_eqs, init_conc, t_grid, k_list, kwargs)

    df = pd.DataFrame({"time": t_grid})
    for i, s in enumerate(parameters["species"]):
        vals_nominal = [nominal(val) for val in sol[:, i]]
        vals_std     = [stdev(val)   for val in sol[:, i]]
        df[s] = vals_nominal
        df[s + "_std"] = vals_std

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

def _safe_spd(P, base=1e-9, tries=8):
    P = torch.nan_to_num(P, nan=0.0, posinf=1e6, neginf=-1e6)
    P = 0.5 * (P + P.transpose(-1, -2))
    n = P.shape[0]
    I = torch.eye(n, dtype=P.dtype, device=P.device)
    jit = base
    for _ in range(tries):
        try:
            return torch.linalg.cholesky(P + jit * I)
        except Exception:
            jit = jit * 10.0
    d = torch.diagonal(P)
    d = torch.nan_to_num(d, nan=jit, posinf=1e6, neginf=jit)
    d = torch.clamp(d, min=jit)
    return torch.linalg.cholesky(torch.diag(d))

def _safe_inv(A, base=1e-9, tries=8):
    A = torch.nan_to_num(A, nan=0.0, posinf=1e6, neginf=-1e6)
    A = 0.5 * (A + A.transpose(-1, -2))
    n = A.shape[0]
    I = torch.eye(n, dtype=A.dtype, device=A.device)
    jit = base
    for _ in range(tries):
        try:
            Ai = torch.linalg.inv(A + jit * I)
            if torch.all(torch.isfinite(Ai)):
                return Ai
        except Exception:
            pass
        jit = jit * 10.0
    d = torch.nan_to_num(torch.diagonal(A), nan=jit, posinf=1e6, neginf=jit)
    d = torch.clamp(d, min=jit)
    return torch.diag(1.0 / d)

def safe_sqrt(x, eps=1e-8):
    return torch.sqrt(torch.clamp(x, min=eps))

_NAN_TRACE_HITS = {}
def _nan_check(name, tsr):
    try:
        bad = (~torch.isfinite(tsr)).any().item()
    except Exception:
        bad = True
    if bad and name not in _NAN_TRACE_HITS:
        _NAN_TRACE_HITS[name] = True
        print(f"  [NaN-TRACE] first non-finite: {name}  shape={tuple(tsr.shape) if hasattr(tsr,'shape') else '?'}")
    return bad

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
        x = torch.nan_to_num(x, nan=0.0, posinf=1e6, neginf=-1e6)   # sanitize state too
        L = _safe_spd(P)
        X = [x]
        scale = torch.sqrt(torch.tensor(n + self.lambda_, dtype=P.dtype, device=P.device))
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
        K = P_xz @ _safe_inv(P_zz)
        x_up = x_pred + K @ (z - z_pred)
        P_up = P_pred - K @ P_zz @ K.T
        P_up = 0.5 * (P_up + P_up.transpose(-1, -2))   # keep symmetric for next step
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
        P_pred = P_pred + Q
        P_pred = 0.5 * (P_pred + P_pred.transpose(-1, -2))
        return x_pred, P_pred
# }}}

# BayesMM{{{
class BayesMM:
    """ The Bayesian Metamodel class.
    Attributes:
      models (list[SurrogateModel]): list of input models
      param_names_list (list[list[str]] | None):
      universal_dt (float): universal timestep
      constraint_maps (list[Callable | None] | None):
          One per model
              f(x_before, x_after, params, dt) -> x_propagated
          where x_before is the state BEFORE the current coupling step,
          x_after is the state AFTER, params is the current parameter dict
          for that model, and dt is the universal_dt.
          The callable function returns a new state tensor in which
          dependent variables have been recomputed.
          Pass None for the whole list to disable the feature entirely.
    """
    def __init__(self, models, param_names_list=None, universal_dt=0.01, constraint_maps=None):
        self.models           = models
        self.param_names_list = param_names_list or [[] for _ in models]
        self.universal_dt     = universal_dt
        self.coupler          = Coupler()
        self.surrogates       = [UnscentedKalmanFilter(model.state_dim, model.Q) for model in models]
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

    def couple_models(self, data_list, connecting_vars):
        couplings = []
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
        else: a = 1.0
        if 'shifts' in c_dict and c_dict['shifts'] is not None and len(c_dict['shifts']) > local_idx:
            b = c_dict['shifts'][local_idx]
        elif 'offsets' in c_dict and c_dict['offsets'] is not None and len(c_dict['offsets']) > local_idx:
            b = c_dict['offsets'][local_idx]
        else: b = 0.0
        return (
            self._as_scalar_tensor(a, device, dtype=dtype),
            self._as_scalar_tensor(b, device, dtype=dtype),
        )

    def _compute_effective_omega(self, sig_t, sig_c, w_j=None):
        sig_t2    = torch.clamp(sig_t**2, min=1e-12)
        sig_c2    = torch.clamp(sig_c**2, min=1e-12)
        omega_opt = sig_t2 / (sig_t2 + sig_c2)
        if getattr(self, "learn_omegas", False) and (w_j is not None):
            omega_eff = torch.clamp(0.5 * w_j + 0.5 * omega_opt, 0.0, 1.0)
        else:
            omega_eff = omega_opt
        return omega_eff

    def _apply_physical_coupling_correction(self, c_dict, local_idx, states, covs, data_list, t):
        """
        Kalman-style residual correction for one participating model:
            y_i  = alpha_i * x_i + beta_i
            r    = g_c({y_i}, t=t)
            K_j  = alpha_j * var_j / (sum_i alpha_i^2 * var_i + sigma_c^2)
            x_j <- x_j - K_j * r
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
        try:
            r = c_dict['g_c'](z_mapped, t=t)
        except TypeError:
            r = c_dict['g_c'](z_mapped)

        r = r.to(device=device, dtype=dtype).reshape(())

        var_name_self = c_dict['vars'][local_idx]
        var_idx_self  = data_list[m_self].columns.get_loc(var_name_self) - 1
        alpha_j       = alphas[local_idx]
        var_j         = vars_all[local_idx]
        K_j           = alpha_j * var_j / total_var

        x_corrected = x_up.clone()
        if not c_dict["Kalman_style"]:
            if local_idx == 0:
                return x_up
            # r is used as a result
            x_corrected[var_idx_self] = r
        else:
            # r is used as a residual
            x_corrected[var_idx_self] = x_corrected[var_idx_self] - K_j * r
        return x_corrected


    def _propagate_state_corrections(self, states_before:list, states_after:list, parameters:list) -> list:
        """
        After coupling steps have corrected one or more primary state
        variables in a model, propagate those corrections to dependent
        variables using each model's constraint_map callable.

        Args:
          states_before(list[Tensor]):
              Per-model state vectors BEFORE any coupling correction this
              timestep (i.e. the raw UKF posterior).
          states_after(list[Tensor]):
              Per-model state vectors AFTER all coupling corrections.  Only
              the entries for models whose constraint_map is not None are
              modified.
          parameters(list[dict]):
              Per-model parameter dicts (same as passed to inference).

        Returns:
          list[Tensor]: Updated states_after with dependent variables recomputed.
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
                x_propagated = cm(x_before, x_after, parameters[i], self.universal_dt)
            except Exception as exc:
                raise RuntimeError(
                    f"constraint_maps[{i}] raised an error during "
                    f"propagation at this timestep: {exc}"
                ) from exc
            result[i] = x_propagated
        return result



    def inference( self, data_list, couplings, steps, params_list=None, M2_interp=None, progress_bar=True, differentiable_params=False, detach_covs=True):
        """
        Multi-model UKF inference with statistical and physical couplings.
        """
        self._coerce_coupling_tensors(couplings)
        self.ensure_maps_on_couplings(couplings)

        n_models   = len(self.models)
        state_dims = [model.state_dim for model in self.models]
        states = [
            torch.tensor(data.iloc[0, 1:1+dim].values, dtype=torch.float32, requires_grad=True)
            for data, dim in zip(data_list, state_dims)
        ]
        covariances = [torch.eye(dim, dtype=torch.float32) * 1.0 for dim in state_dims]

        def _coerce_param_value(v):
            #if isinstance(v, torch.Tensor):
            #    if differentiable_params:
            #        return v if v.requires_grad else v.clone()
            #    return v.detach().clone().to(dtype=torch.float32)

            if isinstance(v, torch.Tensor):
                if differentiable_params:
                    return v.to(dtype=torch.float32)
                return v.detach().clone().to(dtype=torch.float32)

            if isinstance(v, (int, float, np.number)):
                return torch.tensor(float(v), dtype=torch.float32)
            return None


        # Get parameters
        parameters = [{} for _ in self.models]
        if params_list is not None:
            if len(params_list) != n_models:
                raise ValueError(f"params_list length ({len(params_list)}) must equal number of models ({len(self.models)}).")
            parameters = []
            for i, raw_params in enumerate(params_list):
                names = self.param_names_list[i]
                p = {}
                if isinstance(raw_params, dict):
                    for name in names:
                        if name in raw_params:
                            tv = _coerce_param_value(raw_params[name])
                            if tv is not None:
                                p[name] = tv
                else:
                    raise TypeError(f"params_list[{i}] must be a dict, got {type(raw_params)}")
                parameters.append(p)

        # Output storage
        states_over_time      = [[] for _ in self.models]
        param_results         = [[] for _ in self.models]
        covariances_over_time = [[] for _ in self.models]
        z_pred_over_time      = [[] for _ in self.models]
        P_zz_over_time        = [[] for _ in self.models]
        Zi_over_time          = [[] for _ in self.models]
        P_pred_over_time      = [[] for _ in self.models]

        stat_couplings = [c for c in couplings if c.get('coupling_type', 'statistical') != 'physical']
        phys_couplings = [c for c in couplings if c.get('coupling_type') == 'physical']
        coupling_states = [c['mu_c'].clone().detach().requires_grad_(True) for c in stat_couplings]
        coupling_covs = [(c['sigma_c'].clone().detach().requires_grad_(True) ** 2).view(1, 1) for c in stat_couplings]

        progress = tqdm(range(steps), desc="Inference") if progress_bar else range(steps)

        for t in progress:
            _x_pred   = [None] * n_models
            _P_pred   = [None] * n_models
            _z_pred   = [None] * n_models
            _P_zz     = [None] * n_models
            _params_t = [None] * n_models
            new_states = []
            new_covs   = []

            for i, (model, surrogate, state, cov, params) in enumerate(
                zip(self.models, self.surrogates, states, covariances, parameters)
            ):
                t_phys = t * self.universal_dt
                forward_fn = lambda x, dt, _p=params, _t=t_phys: model.f(x, dt, _p, _t)
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
                    obs_std = torch.nan_to_num(obs_std, nan=1e-3, posinf=1e3, neginf=1e-3)
                    obs_std = torch.clamp(obs_std, min=1e-3)
                    R = torch.diag(obs_std ** 2)
                except IndexError:
                    R = model.R(state)

                _nan_check(f"x_pred[m={i},t={t}]", x_pred)
                _nan_check(f"P_pred[m={i},t={t}]", P_pred)
                _nan_check(f"R[m={i},t={t}]", R)
                x_up, P_up, z_pred, P_zz = surrogate.update(x_pred, P_pred, obs_mean, model.g, R)
                _nan_check(f"z_pred[m={i},t={t}]", z_pred)
                _nan_check(f"P_zz[m={i},t={t}]", P_zz)
                _nan_check(f"x_up[m={i},t={t}]", x_up)

                _x_pred[i]   = x_pred
                #_P_pred[i]   = P_pred.detach().clone()
                _P_pred[i] = P_pred.detach().clone() if detach_covs else P_pred
                _z_pred[i]   = z_pred
                _P_zz[i]     = P_zz
                _params_t[i] = {
                    name: val.detach().numpy() if isinstance(val, torch.Tensor) else val
                    for name, val in params.items()
                }

                new_states.append(x_up)
                new_covs.append(P_up)

            states_before_coupling = [s.detach().clone() for s in new_states]

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
                            sig_t2 = torch.clamp(
                                torch.nan_to_num(P_up[var_idx, var_idx], nan=1e-8,
                                                 posinf=1e6, neginf=1e-8), min=1e-8)
                            mu_c   = coupling_states[c_idx]
                            sig_c2 = torch.clamp(
                                torch.nan_to_num(coupling_covs[c_idx].view(-1)[0], nan=1e-8,
                                                 posinf=1e6, neginf=1e-8), min=1e-8)
                            omega_opt = sig_t2 / torch.clamp(sig_t2 + sig_c2, min=1e-6)
                            w_j       = c['weights'][model_local]
                            if getattr(self, "learn_omegas", False):
                                omega_eff = torch.clamp(0.5 * w_j + 0.5 * omega_opt, 0.0, 1.0)
                            else:
                                omega_eff = omega_opt

                            x_blend        = omega_eff * mu_c + (1.0 - omega_eff) * mu_t
                            x_up2          = x_up.clone()
                            x_up2[var_idx] = x_blend
                            x_up           = x_up2

                    new_states[i] = x_up

            if phys_couplings:
                base_states = [s.clone() for s in new_states]
                cur_covs    = new_covs
                accum       = [torch.zeros_like(s) for s in new_states]
                for c in phys_couplings:
                    for local_idx, m in enumerate(c['models']):
                        corrected = self._apply_physical_coupling_correction(
                            c_dict    = c,
                            local_idx = local_idx,
                            states    = base_states,   # always the snapshot
                            covs      = cur_covs,
                            data_list = data_list,
                            t         = t,
                        )
                        accum[m] = accum[m] + (corrected - base_states[m])
                merged = []
                for m in range(n_models):
                    s = base_states[m] + accum[m]
                    # Guard: a single non-finite correction must not poison the
                    # rest of the trajectory or the downstream MultivariateNormal.
                    s = torch.nan_to_num(s, nan=0.0, posinf=1e6, neginf=-1e6)
                    merged.append(s)
                new_states = merged


            if any(cm is not None for cm in self.constraint_maps):
#                print("new_states(before) = ",new_states)
                new_states = self._propagate_state_corrections(
                    states_before = states_before_coupling,
                    states_after  = new_states,
                    parameters    = parameters,
                )
#                print("new_states(after) = ",new_states)

            for i in range(n_models):
                states_over_time[i].append(new_states[i])
                covariances_over_time[i].append(new_covs[i].detach().clone() if detach_covs else new_covs[i])
                param_results[i].append(_params_t[i])
                Zi_over_time[i].append(_x_pred[i])
                P_pred_over_time[i].append(_P_pred[i])
                z_pred_over_time[i].append(_z_pred[i])
                P_zz_over_time[i].append(_P_zz[i])

            states      = new_states
            covariances = new_covs

            # Recompute statistical priors for next timestep
            postfix_info = {}
            for c_idx, c in enumerate(stat_couplings):
                mus, sigmas = [], []
                for model_local, m in enumerate(c['models']):
                    var_name = c['vars'][model_local]
                    var_idx  = data_list[m].columns.get_loc(var_name) - 1
                    mu       = states[m][var_idx]
                    sigma    = safe_sqrt(covariances[m][var_idx, var_idx], eps=1e-8)
                    mus.append(mu)
                    sigmas.append(sigma)
                mu_c, sigma_c               = self.coupler.mixture_prior(mus, sigmas)
                coupling_states[c_idx]      = mu_c
                coupling_covs[c_idx]        = (sigma_c ** 2).view(1, 1)
                postfix_info[f"cov_c{c_idx}"] = f"{sigma_c**2:.4e}"

            if progress_bar and postfix_info:
                progress.set_postfix(postfix_info)
        states_over_time      = [torch.stack(s) for s in states_over_time]
        P_pred_over_time      = [torch.stack(p) for p in P_pred_over_time]
        covariances_over_time = [torch.stack(c) for c in covariances_over_time]
        return (states_over_time, param_results, covariances_over_time,
                z_pred_over_time, P_zz_over_time, Zi_over_time,P_pred_over_time)


    def compute_data_likelihood_energy(self, steps, z_pred_over_time, P_zz_over_time, data_list, delta=3.0):
        total = torch.zeros((), dtype=torch.float64)
        for i, (z_pred_list, P_zz_list, data_df) in enumerate(
                zip(z_pred_over_time, P_zz_over_time, data_list)):
            state_dim = self.models[i].state_dim
            obs_names = data_df.columns[1:1+state_dim]
            obs_mean  = torch.tensor(np.stack([data_df[n] for n in obs_names], axis=1),
                                     dtype=torch.float64)
            for t in range(steps):
                loc = torch.nan_to_num(z_pred_list[t].to(torch.float64),
                                       nan=0.0, posinf=1e6, neginf=-1e6)
                cov = torch.nan_to_num(P_zz_list[t].to(torch.float64),
                                       nan=1.0, posinf=1e6, neginf=1e-6)
                v   = torch.clamp(torch.diagonal(cov), min=1e-12)
                z   = (obs_mean[t] - loc) / torch.sqrt(v)
                az  = torch.abs(z)
                quad = torch.where(az <= delta, 0.5 * z * z,
                                   delta * az - 0.5 * delta * delta)
                total = total + torch.sum(quad + 0.5 * torch.log(v))
        return total


    def compute_state_variable_likelihood_energy(self, steps, states_over_time, Zi_over_time, P_pred_over_time):
        total_nll_states = 0.0
        for x_up_list, x_pred_list, P_pred_list in zip(states_over_time, Zi_over_time, P_pred_over_time):
            for t in range(steps):
                _loc = torch.nan_to_num(x_pred_list[t], nan=0.0, posinf=1e6, neginf=-1e6)
                _cov = torch.nan_to_num(P_pred_list[t], nan=0.0, posinf=1e6, neginf=-1e6)
                _cov = 0.5 * (_cov + _cov.transpose(-1, -2))
                _cov = _cov + 1e-6 * torch.eye(_cov.shape[0], dtype=_cov.dtype, device=_cov.device)
                dist = MultivariateNormal(_loc, covariance_matrix=_cov)
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
                raise ValueError(f"Unexpected tensor shape for Zi_over_time[{m}]: {z_m.shape}")
        else:
            z_t = z_m[t]
            if isinstance(z_t, torch.Tensor):
                return z_t[var_idx].reshape(()).to(device=device, dtype=dtype)
            return self._as_scalar_tensor(z_t[var_idx], device, dtype=dtype)

    def _get_cov_scalar(self, state_covs_over_time, m, t, var_idx, device,dtype=torch.float64, fallback_var=1e-3):
        if state_covs_over_time is None:
            return torch.tensor(fallback_var, device=device, dtype=dtype).reshape(())
        P_m = state_covs_over_time[m]
        if isinstance(P_m, torch.Tensor):
            if P_m.dim() == 3:
                return torch.clamp(P_m[t, var_idx, var_idx], min=1e-12).reshape(()).to(device=device, dtype=dtype)
            raise ValueError(f"Unexpected tensor shape for state_covs_over_time[{m}]: {P_m.shape}")
        P_t = P_m[t]
        if isinstance(P_t, torch.Tensor):
            return torch.clamp(P_t[var_idx, var_idx], min=1e-12).reshape(()).to(device=device, dtype=dtype)
        val = P_t[var_idx][var_idx]
        return torch.clamp(self._as_scalar_tensor(val, device, dtype=dtype), min=1e-12)

    def _apply_coupling_map_and_slack(self, c_dict, mu_j, var_j, local_idx, tau2_c):
        alpha, beta = self._get_map_params(c_dict, local_idx, mu_j.device, mu_j.dtype)
        y_mu        = alpha * mu_j + beta
        y_var       = torch.clamp(alpha * alpha * var_j + tau2_c, min=1e-12)
        return y_mu, y_var


    def _gather_mu_var_for_coupling(self, t, c_dict, states_over_time,
        state_covs_over_time, data_list, device, tau2_c, fallback_var=1e-3):
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


    def _gaussian_energy(self, d, veff, delta=3.0):
        """
        Bounded per-term Gaussian energy: 0.5*z^2 + 0.5*log(veff), with the
        quadratic replaced by z = d/sqrt(veff).
        """
        veff = torch.clamp(veff.reshape(-1), min=1e-12)
        z    = d.reshape(-1) / torch.sqrt(veff)
        az   = torch.abs(z)
        quad = torch.where(az <= delta, 0.5 * z * z, delta * az - 0.5 * delta * delta)
        logdet = 0.5 * torch.log(veff)
        return torch.sum(quad + logdet)

    def _coupling_difference_and_var(self, c, cur_states, cur_covs, data_list, device, dtype, t):
        is_phys   = c.get('coupling_type', 'statistical') == 'physical'
        is_value  = is_phys and not c.get('Kalman_style', True)

        y_mus, y_vars = [], []
        for j, m in enumerate(c['models']):
            var_idx = data_list[m].columns.get_loc(c['vars'][j]) - 1
            mu_j  = cur_states[m][var_idx].to(device=device, dtype=dtype)
            var_j = torch.clamp(cur_covs[m][var_idx, var_idx].to(device, dtype), min=1e-12)
            a, b  = self._get_map_params(c, j, device, dtype)
            y_mus.append(a * mu_j + b)
            y_vars.append(torch.clamp(a * a * var_j, min=1e-12))
        y = torch.stack(y_mus)   # (J,)
        v = torch.stack(y_vars)  # (J,)

        # 0 for statistical; sigma_c^2 for physical
        if is_phys:
            sigma_c = self._as_scalar_tensor(c['sigma_c'], device, dtype)
            tau2    = torch.clamp(sigma_c ** 2, min=0.0)
        else:
            tau2    = torch.zeros((), device=device, dtype=dtype)

        if is_value:
            try:    g = c['g_c'](list(y_mus), t=t)
            except TypeError: g = c['g_c'](list(y_mus))
            g = g.to(device=device, dtype=dtype).reshape(())
            # The discrepancy y_j - g carries the uncertainty of BOTH the target
            # and the source
            v_total = torch.sum(v) + tau2
            d = (y - g).reshape(-1)
            # only the TARGET participants carry a residual (source d = 0)
            return d, v_total.expand_as(d)

        # statistical AND physical-difference: identical consensus construction.
        # tau2 is the ONLY difference between them, and it is 0 for statistical.
        prec    = 1.0 / v
        mu_pool = torch.sum(prec * y) / torch.sum(prec)
        d       = y - mu_pool
        veff    = v + 1.0 / torch.sum(prec) + tau2
        return d, veff


    def compute_loss(self, data_list, couplings, steps, Zi_over_time, states_over_time,
                     covariances_over_time, z_pred_over_time, P_zz_over_time,
                     P_pred_over_time, device=torch.device("cpu"), dtype=torch.float64,
                     include_states=False, reduce_mean=False, delta=3.0):
        """
        Harmony loss, differentiable, returns a torch scalar
        """
        self._coerce_coupling_tensors(couplings)
        self.ensure_maps_on_couplings(couplings)
        n_models       = len(self.models)
        n_couplings    = len(couplings)
        total_obs_data = sum(steps * self.models[i].state_dim for i in range(n_models))

        nll_cpl = torch.zeros((), dtype=dtype, device=device)
        for t in range(steps):
            cur_states = [states_over_time[m][t] for m in range(n_models)]
            cur_covs   = [covariances_over_time[m][t] for m in range(n_models)]
            for c in couplings:
                d, veff = self._coupling_difference_and_var(
                    c, cur_states, cur_covs, data_list, device, dtype, t)
                nll_cpl = nll_cpl + self._gaussian_energy(d, veff)

        nll_data   = self.compute_data_likelihood_energy(
            #steps, z_pred_over_time, P_zz_over_time, data_list, delta=2.0).to(dtype)
            steps, z_pred_over_time, P_zz_over_time, data_list, delta=delta).to(dtype)
        nll_states = self.compute_state_variable_likelihood_energy(
            steps, states_over_time, Zi_over_time, P_pred_over_time).to(dtype) \
            if include_states else torch.zeros((), dtype=dtype, device=device)

        total = nll_cpl + nll_data #+ nll_states
        if reduce_mean:
            total = total / total_obs_data

        # cache normalized pieces for logging/back-compat
        self.num_obs      = total_obs_data
        self.n_couplings  = n_couplings
        self.nll_cpl_n    = float(nll_cpl / total_obs_data)
        self.nll_data_n   = float(nll_data / total_obs_data)
        self.nll_states_n = float(nll_states / total_obs_data)
        self.nll_cpl      = nll_cpl
        self.nll_data     = nll_data
        self.nll_states   = nll_states
        self.nll          = total
        return total


# }}}

# Physical coupling:{{{

@dataclass
class PhysicalCoupling:
    """Declares a single physical coupling constraint and serializes to a dict."""
    model_indices : List[int]
    var_names     : List[str]
    g_c           : Callable
    sigma_c       : float = 0.1
    lambdas       : List[float] = field(default_factory=list)
    shifts        : List[float] = field(default_factory=list)
    weights       : List       = field(default_factory=list)
    Kalman_style  : bool = True

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
            'Kalman_style': self.Kalman_style,
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
    Kalman_style: bool = True
) -> dict:
    """
    Build a physical coupling dict compatible with BayesMM.inference().
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
        Kalman_style  = Kalman_style
    )
    return pc.as_dict()
# }}}

# Helper functions:{{{
def _ordered_values(model: dict, field: str):
    names = list(model.get("param_names", []))
    vals = model.get(field, None)
    if vals is None:
        raise KeyError(f"model['{field}'] is missing")

    if isinstance(vals, dict):
        missing = [n for n in names if n not in vals]
        if missing:
            raise KeyError(f"Missing {field} values for: {missing}")
        return [vals[n] for n in names]

    if isinstance(vals, (list, tuple, np.ndarray)):
        if len(vals) != len(names):
            raise ValueError(
                f"param_names has {len(names)} entries but {field} has {len(vals)} values."
            )
        return list(vals)

    raise TypeError(f"Unsupported type for model['{field}']: {type(vals)}")

def extract_parameter_objects(models_data, requires_grad=False, dtype=torch.float32, device=None):
    params_list, param_names_list = [], []
    for m in models_data:
        names, vals = m["params"]["param_names"], m["params"]["params"]
        params_list.append({
            n: torch.tensor(float(v), dtype=dtype, device=device, requires_grad=requires_grad)
            for n, v in zip(names, vals)
        })
        param_names_list.append(names)
    return params_list, param_names_list

def extract_parameter_std_objects(models_data):
    return [
        {n: (None if s is None else float(s)) for n, s in zip(m["params"]["param_names"], m["params"]["params_std"])}
        for m in models_data
    ]


def _torch_to_numpy_state(x_torch: torch.Tensor) -> np.ndarray:
    return x_torch.detach().cpu().numpy().astype(float)

def _torch_params_to_numpy_list(params, param_names):
    return [float(params[name].detach().cpu().item()) for name in param_names]

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


def generate_surrogate(observables, model_parameters, ode_fn, sim_time, dt,
    initial_noise_scale=0.01, transition_cov_scale=0.01, emission_cov_scale=1.0,
    observation_noise_scale=0.001, noise_model_type='time-invariant',ode_kwargs: dict = None):
    """
    """
    obs_mean  = observables[:, :, 0]
    obs_std   = observables[:, :, 1]
    steps     = obs_mean.shape[0]
    state_dim = len(model_parameters["InitConc"])

    param_names = list(model_parameters["param_names"])
    param_vals  = _ordered_values(model_parameters, "params")
    param_dict  = {
        name: torch.tensor(val, dtype=torch.float32, requires_grad=True)
        for name, val in zip(param_names, param_vals)
    }

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
                x * torch.normal(mean=emission_cov_scale, std=emission_cov_scale * 1e-2, size=x.shape)
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

    bayes_mm = BayesMM(models=[surrogate], param_names_list=[param_names], universal_dt=dt)
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


def _powerset(iterable):
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def auto_identify_couplings(data_dfs, params_list, models_data, sim_time,
                            params_list_for_inf=None, candidate_phys_couplings=None):
    """
    Joint search over statistical coupling candidates & physical coupling subsets.
    Returns the configuration (stat_vars, phys_subset, Zi_over_time) that
    minimizes the normalized total loss.
    """
    candidate_phys_couplings = candidate_phys_couplings or []
    phys_subsets             = _powerset(candidate_phys_couplings)

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
        bayes_mm._coerce_coupling_tensors(all_couplings)
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

        loss = bayes_mm.compute_loss(data_dfs, all_couplings, steps,
            Zi_over_time, states_over_time, covariances_over_time,
            z_pred_over_time, P_zz_over_time, P_pred_over_time, device=device, dtype=dtype)

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



def build_metamodel(
    models_data,
    sim_time,
    connecting_vars          = None,
    couplings                = None,
    candidate_phys_couplings = None,
    candidate_graphs         = None,
    run_param_opt            = True,
    lr                       = 4e-2,
    max_iterations           = 50,
    plot_input_models        = False
):
    """
    Construct the Bayesian metamodel.

    Args:
      models_data(list[dict]):
          One dict per model.  Required keys: 'obs_model', 'surrogate', 'params',
          'ode_fn', 'dt', 'transition_cov_scale', 'emission_cov_scale',
          'observation_noise_scale', 'noise_model_type', 'param_names'.
      sim_time(float): Total simulation time.
      connecting_vars(list | None): Statistical connecting variable pairs.  None -> auto-search.
      couplings(list | None): Pre-built statistical coupling dicts.  None -> built from connecting_vars.
      candidate_phys_couplings(list[dict] | None): Physical coupling dicts from make_physical_coupling().
    """
    obs_models       = [m['obs_model']   for m in models_data]
    params_list      = [m['params']      for m in models_data]
    ode_fns          = [m['ode_fn']      for m in models_data]
    dts              = [m['dt']          for m in models_data]
    param_names_list = [m['param_names'] for m in models_data]
    universal_dt   = min(dts)
    steps          = int(sim_time / universal_dt)
    universal_time = np.arange(0, sim_time, universal_dt)[:steps]

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

    # Build surrogates
    for idx, (obs_processed, params, ode_fn) in enumerate(zip(processed_obs, params_list, ode_fns)):
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

    #Coupling selection
    selected_phys_couplings = candidate_phys_couplings or []

    if connecting_vars is None:
        params_list_for_inf = [
            {name: params['params'][i] for i, name in enumerate(names)}
            for params, names in zip(params_list, param_names_list)
        ]
        connecting_vars, selected_phys_couplings, _ = auto_identify_couplings(
            data_dfs, params_list, models_data, sim_time,
            params_list_for_inf      = params_list_for_inf,
            candidate_phys_couplings = candidate_phys_couplings or [],
            candidate_graphs         = candidate_graphs,
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
    bayes_mm.all_couplings = all_couplings
    bayes_mm._coerce_coupling_tensors(all_couplings)
    bayes_mm.ensure_maps_on_couplings(all_couplings)

    # Final inference
    params_list_for_opt = [
        {name: params['params'][i] for i, name in enumerate(names)}
        for params, names in zip(params_list, param_names_list)
    ]

    (states_over_time, param_results, covariances_over_time, z_pred_over_time,
     P_zz_over_time, Zi_over_time, P_pred_over_time) = bayes_mm.inference(
        data_list   = data_dfs,
        couplings   = all_couplings,
        steps       = steps,
        params_list = params_list_for_opt,
    )

    bayes_mm.loss = loss = bayes_mm.compute_loss(
        data_dfs, all_couplings, steps,
        Zi_over_time, states_over_time, covariances_over_time,
        z_pred_over_time, P_zz_over_time, P_pred_over_time,
        reduce_mean=True
    ).float()
    bayes_mm.loss_std = 0.0
    print(f"BMM loss       = {loss:.4f}")
    print(f"  nll_cpl/obs  = {bayes_mm.nll_cpl_n:.4f}")
    print(f"  nll_data/obs = {bayes_mm.nll_data_n:.4f}")

    states = [np.array(r.detach().numpy()) for r in states_over_time]

    coupled_dfs = []
    for idx, (state, params) in enumerate(zip(states, params_list)):
        coupled_dfs.append(pd.DataFrame({
            'time': universal_time,
            **{name: state[:, i] for i, name in enumerate(params["species"])},
        }))

    figs = []
    if plot_input_models:
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

# laplace_inference:{{{

def make_torch_ode_forward(ode_fn: Callable, param_names: Sequence[str],
        ode_kwargs: Optional[dict] = None, t_arg: float = 0.0) -> Callable:
    """
    Differentiable replacement for the numpy `forward_fn` bridge.

    Parameters
    ----------
    ode_fn : Callable
        Signature ode_fn(x, t, k_list, **ode_kwargs) -> sequence of dx/dt.
    param_names : Sequence[str]
        Ordering used to build k_list from the params dict.
    ode_kwargs : dict, optional
        Extra keyword arguments forwarded to ode_fn.
    t_arg : float
        Value passed as the `t` argument of ode_fn.

    Returns
    -------
    Callable with signature forward(x, dt_step, params, t=None) matching
    SurrogateModel.f as required inside BayesMM.inference.
    """
    ode_kwargs = dict(ode_kwargs or {})
    names = list(param_names)

    def forward(x: torch.Tensor, dt_step: float, params: dict, t=None) -> torch.Tensor:
        k_list = [params[n] for n in names]
        dxdt = ode_fn(x, t_arg, k_list, **ode_kwargs)
        if isinstance(dxdt, torch.Tensor):
            d = dxdt.reshape(x.shape).to(dtype=x.dtype, device=x.device)
        else:
            elems = []
            for e in dxdt:
                if isinstance(e, torch.Tensor):
                    elems.append(e.reshape(()).to(dtype=x.dtype, device=x.device))
                else:
                    val = e.nominal_value if hasattr(e, "nominal_value") else float(e)
                    elems.append(torch.as_tensor(float(val), dtype=x.dtype,
                                                 device=x.device))
            d = torch.stack(elems).reshape(x.shape)
        return x + dt_step * d

    return forward


def torchify_forwards(bayes_mm, models_data: List[dict], model_indices: Optional[Sequence[int]] = None) -> List[int]:
    """
    Replace `bayes_mm.models[i].f` with a torch-native forward for each model.
    Models with ode_fn=None and no 'torch_forward' are left alone.
    Returns the list of indices that were modified.
    """
    touched = []
    for i, m in enumerate(models_data):
        if model_indices is not None and i not in model_indices:
            continue
        override = m.get("torch_forward", None)
        if override is not None:
            bayes_mm.models[i].f = override
            touched.append(i)
            continue
        if m.get("ode_fn", None) is None:
            continue
        bayes_mm.models[i].f = make_torch_ode_forward(
            m["ode_fn"], m["param_names"], m.get("ode_kwargs", None)
        )
        touched.append(i)
    return touched

def _inference_kwargs(bayes_mm) -> dict:
    try:
        sig = inspect.signature(bayes_mm.inference)
    except (TypeError, ValueError):
        return {}
    return {"detach_covs": False} if "detach_covs" in sig.parameters else {}

def _run_inference(bayes_mm, data_list, couplings, steps, plist, extra_kw):
    return bayes_mm.inference(
        data_list=data_list,
        couplings=couplings,
        steps=steps,
        params_list=plist,
        progress_bar=False,
        differentiable_params=True,
        **extra_kw, )

def _metamodel_nll(bayes_mm, data_list, couplings, steps, out, device, dtype,
                   objective: str = "marginal"):
    """
    -log p(D | theta) from the filter.  Used by objective="marginal"/"harmony".
    "marginal": the UKF prediction error decomposition, un-normalized.
    "harmony": compute_loss_diff * total_obs.  The harmony loss is a coupling graph selection score.
    """
    (states_ot, _, cov_ot, zpred_ot, Pzz_ot, Zi_ot, Ppred_ot) = out

    if objective == "marginal":
        return bayes_mm.compute_data_likelihood_energy(
            steps, zpred_ot, Pzz_ot, data_list).to(dtype)

    if objective == "harmony":
        total_obs = float(sum(steps * m.state_dim for m in bayes_mm.models))
        return bayes_mm.compute_loss(
            data_list, couplings, steps, Zi_ot, states_ot, cov_ot,
            zpred_ot, Pzz_ot, Ppred_ot, device=device, dtype=dtype).to(dtype) #* total_obs

    raise ValueError(f"objective must be trajectory|marginal|harmony, got {objective!r}")


def harmony_loss_at_posterior(bayes_mm, data_list, couplings, steps, params_list,
    theta_map, cov, free_keys, reduce_mean=True, device=torch.device("cpu"), dtype=torch.float64):
    """
    Harmony loss and its uncertainty, propagated from the parameter posterior.
    Returns (L_hat, sigma_L).
    """
    bayes_mm._coerce_coupling_tensors(couplings)
    bayes_mm.ensure_maps_on_couplings(couplings)

    fixed = [
        {n: (v.detach().clone().to(device=device, dtype=dtype) if torch.is_tensor(v)
             else torch.tensor(float(v), dtype=dtype, device=device))
         for n, v in p.items()}
        for p in params_list
    ]

    def L_of_theta(tv):
        plist = [dict(f) for f in fixed]
        for i, (m, name) in enumerate(free_keys):
            plist[m][name] = tv[i]
        out = bayes_mm.inference(
            data_list=data_list, couplings=couplings, steps=steps,
            params_list=plist, progress_bar=False,
            differentiable_params=True, detach_covs=False)
        (states_ot, _, cov_ot, zpred_ot, Pzz_ot, Zi_ot, Ppred_ot) = out
        return bayes_mm.compute_loss(
            data_list, couplings, steps, Zi_ot, states_ot, cov_ot,
            zpred_ot, Pzz_ot, Ppred_ot, device=device, dtype=dtype,
            reduce_mean=reduce_mean)

    tv = torch.tensor([theta_map[k] for k in free_keys],
                      dtype=dtype, device=device, requires_grad=True)
    L = L_of_theta(tv)
    gL, = torch.autograd.grad(L, tv, allow_unused=True)
    if gL is None:
        return float(L.detach()), 0.0     # loss doesn't depend on phi -> no propagated uncertainty
    gL = torch.nan_to_num(gL.detach(), nan=0.0)
    var_L = float(gL @ cov.to(gL.dtype) @ gL)
    return float(L.detach()), float(np.sqrt(max(var_L, 0.0)))


def build_trajectory_targets(bayes_mm, data_list, couplings, steps, params_list):
    """One inference run at params_list.  Returns (states, covs), detached."""
    plist = [
        {n: (v.detach().clone() if torch.is_tensor(v)
             else torch.tensor(float(v))).to(torch.float32)
         for n, v in p.items()}
        for p in params_list
    ]
    with torch.no_grad():
        out = bayes_mm.inference(data_list=data_list, couplings=couplings,
                                 steps=steps, params_list=plist, progress_bar=False)
    return [s.detach() for s in out[0]], [c.detach() for c in out[2]]


def _trajectory_nll(bayes_mm, targets, weights, plist, free_models, dt, dtype,
                    profile_scale=True):
    total = None
    for m in free_models:
        X, W = targets[m], weights[m]
        f = bayes_mm.models[m].f
        r = torch.stack([X[t + 1] - f(X[t], dt, plist[m], t * dt) for t in range(X.shape[0] - 1)])
        wr2 = W * r * r
        if profile_scale:
            n = wr2.numel()
            term = 0.5 * n * torch.log(torch.clamp(wr2.mean(), min=1e-300))
        else:
            term = 0.5 * torch.sum(wr2)
        total = term if total is None else total + term
    if total is None:
        raise ValueError("No model with free parameters; nothing to fit.")
    return total.to(dtype)


def _trajectory_chi2(bayes_mm, targets, weights, plist, free_models, dt, n_free, profile_scale=True):
    """
    Weighted chi^2 per dof at theta.
    """
    num, n = 0.0, 0
    scales = []
    with torch.no_grad():
        for m in free_models:
            X, W = targets[m], weights[m]
            f = bayes_mm.models[m].f
            r = torch.stack([X[t + 1] - f(X[t], dt, plist[m], t * dt) for t in range(X.shape[0] - 1)])
            wr2 = W * r * r
            s2 = float(wr2.mean()) if profile_scale else 1.0
            scales.append(s2)
            num += float(wr2.sum()) / max(s2, 1e-300)
            n += wr2.numel()
    return num / max(n - n_free, 1), scales


def _posterior_cov_from_hessian(H: torch.Tensor, prior_prec: torch.Tensor):
    """
    H = H_likelihood + diag(prior_prec).
    """
    H = 0.5 * (H + H.T)
    P = torch.diag(prior_prec.to(H.dtype))
    H_lik = 0.5 * ((H - P) + (H - P).T)
    w, V = torch.linalg.eigh(H_lik)
    n_neg = int((w < 0).sum())
    H_reg = V @ torch.diag(torch.clamp(w, min=0.0)) @ V.T + P
    cov = torch.linalg.inv(0.5 * (H_reg + H_reg.T))
    return 0.5 * (cov + cov.T), w, n_neg


def laplace_parameter_inference(bayes_mm, data_list, couplings,steps,
    params_list,                 # list[dict[name -> value]]: the MAP seed
    param_std_list=None,         # priors; std None/<=0 means "fixed"
    extra_evidence_fn=None,      # e.g. the delta-G Student-t term
    objective: str = "trajectory",   # "trajectory" | "marginal" | "harmony"
    targets=None,                # list[(T,dim)] coupled states; None -> one seed run
    target_covs=None,            # list[(T,dim,dim)] coupled covariances
    traj_scale: float = 1.0,     # multiplier on S; see chi2/dof in the printout
    traj_profile_scale: bool = True,   # profile out the per-model noise scale
    traj_var_floor: float = 1e-10, lr: float = 2e-2, n_iters: int = 400,
    prior_strength: float = 1.0,patience: int = 40,
    grad_rtol: float = 1e-4,     # converged when |g| < grad_rtol * |g_seed|
    grad_tol: Optional[float] = None, rel_tol: float = 1e-9,
    refine_lbfgs: bool = True, lbfgs_iters: int = 200,
    device=torch.device("cpu"), dtype=torch.float64,
    grad_mode: str = "hybrid",   # "autograd" | "fd" | "hybrid"
    hessian_mode: str = "fd",    # "fd" | "autograd"
    fd_rel: float = 1e-3, fd_abs: float = 1e-8,
    fd_rel_hess: Optional[float] = None, verbose: bool = True,):
    """
    MAP + Laplace posterior over the free parameters of every model.

        \mathcal{L}(theta) = nll(theta) + 0.5*sum_i prior_prec_i (theta_i - mu_i)^2
                 + extra_evidence_fn(params, inference_output)

    objective="trajectory" (default)
        Fit theta so each model's dynamics reproduce the coupled posterior trajectory
    objective="marginal" / "harmony"
        See _metamodel_nll.

    Alternate with refit_at_map: fit theta to the trajectory, refit the
    trajectory at theta, repeat.  Pass the refit states back in as `targets`.

    Returns theta_map, theta_std, cov, free_keys, info.
    """
    if objective not in ("trajectory", "marginal", "harmony"):
        raise ValueError(f"objective must be trajectory|marginal|harmony, got {objective!r}")
    if grad_mode not in ("autograd", "fd", "hybrid"):
        raise ValueError(f"grad_mode must be autograd|fd|hybrid, got {grad_mode!r}")
    if hessian_mode not in ("fd", "autograd"):
        raise ValueError(f"hessian_mode must be fd|autograd, got {hessian_mode!r}")

    fd_rel_hess = fd_rel_hess if fd_rel_hess is not None else max(fd_rel, 5e-3)

    bayes_mm._coerce_coupling_tensors(couplings)
    bayes_mm.ensure_maps_on_couplings(couplings)
    extra_kw = _inference_kwargs(bayes_mm)
    dt = float(bayes_mm.universal_dt)

    # free parameters
    free_keys: List[Tuple[int, str]] = []
    theta0, prior_mu_l, prior_prec_l = [], [], []
    for m, params in enumerate(params_list):
        stds = (param_std_list[m] if param_std_list is not None else {}) or {}
        for name, val in params.items():
            v = float(val.detach()) if torch.is_tensor(val) else float(val)
            s = stds.get(name, None)
            if s is None or float(s) <= 0.0:
                continue
            free_keys.append((m, name))
            theta0.append(v)
            prior_mu_l.append(v)
            prior_prec_l.append(prior_strength / float(s) ** 2)

    if not free_keys:
        raise ValueError("No free parameters (all priors are None/<=0).")

    D = len(free_keys)
    free_models = sorted({m for m, _ in free_keys})
    theta = torch.tensor(theta0, dtype=dtype, device=device, requires_grad=True)
    prior_mu = torch.tensor(prior_mu_l, dtype=dtype, device=device)
    prior_prec = torch.tensor(prior_prec_l, dtype=dtype, device=device)

    fixed_snapshot = [
        {name: (val.detach().clone().to(device=device, dtype=dtype)
                if torch.is_tensor(val)
                else torch.tensor(float(val), dtype=dtype, device=device))
         for name, val in params.items()}
        for params in params_list
    ]

    def assemble_params(theta_vec: torch.Tensor) -> List[dict]:
        out = [dict(snap) for snap in fixed_snapshot]
        for idx, (m, name) in enumerate(free_keys):
            out[m][name] = theta_vec[idx]
        return out

    T_list = W_list = None
    if objective == "trajectory":
        if targets is None or target_covs is None:
            if verbose:
                print("[laplace] no targets given; one inference run at the seed")
            targets, target_covs = build_trajectory_targets(
                bayes_mm, data_list, couplings, steps, params_list)
        T_list, W_list = [], []
        for X, P in zip(targets, target_covs):
            Xd = X.detach().to(device=device, dtype=dtype)
            d = torch.diagonal(P.detach().to(device=device, dtype=dtype),
                               dim1=-2, dim2=-1)                  # (T, dim)
            S = float(traj_scale) * (d[:-1] + d[1:])
            T_list.append(Xd)
            W_list.append(1.0 / torch.clamp(S, min=traj_var_floor))

    # the objective
    def _nlp_tensor(theta_vec: torch.Tensor) -> torch.Tensor:
        plist = assemble_params(theta_vec)
        if objective == "trajectory":
            nll, out = _trajectory_nll(bayes_mm, T_list, W_list, plist,
                                       free_models, dt, dtype,
                                       profile_scale=traj_profile_scale), None
        else:
            out = _run_inference(bayes_mm, data_list, couplings, steps, plist, extra_kw)
            nll = _metamodel_nll(bayes_mm, data_list, couplings, steps, out,
                                 device, dtype, objective=objective)
        total = nll.to(dtype) + 0.5 * torch.sum(prior_prec * (theta_vec - prior_mu) ** 2)
        if extra_evidence_fn is not None:
            e = extra_evidence_fn(plist, out)
            if not torch.is_tensor(e):
                e = torch.tensor(float(e), dtype=dtype, device=device)
            # Force scalar
            e = e.to(dtype).sum() if e.dim() > 0 else e.to(dtype)
            total = total + e
        return total.reshape(())

    def _nlp_value(theta_np) -> float:
        with torch.no_grad():
            tv = torch.as_tensor(np.asarray(theta_np, dtype=np.float64),
                                 dtype=dtype, device=device)
            return float(_nlp_tensor(tv))

    def _autograd_grad(theta_vec: torch.Tensor):
        tv = theta_vec.detach().clone().requires_grad_(True)
        loss = _nlp_tensor(tv)
        g, = torch.autograd.grad(loss, tv, allow_unused=True)
        if g is None:
            g = torch.zeros_like(tv)
        return float(loss.detach()), torch.nan_to_num(g.detach(), nan=0.0)

    def _fd_step(x_i: float, rel: float) -> float:
        return rel * max(abs(x_i), fd_abs)

    def _fd_partial(theta_np, i: int, rel: float) -> float:
        h = _fd_step(theta_np[i], rel)
        tp = theta_np.copy(); tp[i] += h
        tm = theta_np.copy(); tm[i] -= h
        return (_nlp_value(tp) - _nlp_value(tm)) / (2.0 * h)

    sens: Dict[Tuple[int, str], float] = {}
    if grad_mode == "fd":
        fd_idx = list(range(D))
    elif grad_mode == "autograd":
        fd_idx = []
    else:
        if objective == "trajectory":
            probes, plist_p = [], []
            for m, params in enumerate(params_list):
                d = {}
                for name, val in params.items():
                    v = float(val.detach()) if torch.is_tensor(val) else float(val)
                    tp = torch.tensor(v, dtype=dtype, device=device, requires_grad=True)
                    d[name] = tp
                    probes.append(((m, name), tp))
                plist_p.append(d)
            nll_p = _trajectory_nll(bayes_mm, T_list, W_list, plist_p,
                                    free_models, dt, dtype,
                                    profile_scale=traj_profile_scale)
            gs = torch.autograd.grad(nll_p, [t for _, t in probes], allow_unused=True)
            sens = {k: (0.0 if g is None else float(torch.nan_to_num(g.abs())))
                    for (k, _), g in zip(probes, gs)}
        else:
            sens = param_sensitivity(bayes_mm, data_list, couplings, steps,
                                     params_list, device=device, dtype=dtype)
        fd_idx = [i for i, k in enumerate(free_keys) if sens.get(k, 0.0) == 0.0]

    if verbose:
        print(f"[laplace] objective={objective} | {D} free parameters: {free_keys}")
        for k in free_keys:
            if k in sens:
                tag = " -> finite differences" if sens[k] == 0.0 else ""
                print(f"[laplace]   |d(nll)/d{k}| = {sens[k]:.3e}{tag}")
        print(f"[laplace] fd coords: {[free_keys[i] for i in fd_idx] or 'none'}")

    def _grad(theta_vec: torch.Tensor):
        if len(fd_idx) == D:
            tnp = theta_vec.detach().cpu().numpy().astype(np.float64)
            base = _nlp_value(tnp)
            g = np.zeros(D, dtype=np.float64)
            for i in range(D):
                g[i] = _fd_partial(tnp, i, fd_rel)
            return base, torch.tensor(g, dtype=dtype, device=device)
        base, g = _autograd_grad(theta_vec)
        if fd_idx:
            tnp = theta_vec.detach().cpu().numpy().astype(np.float64)
            g = g.clone()
            for i in fd_idx:
                g[i] = _fd_partial(tnp, i, fd_rel)
        return base, g

    seed_loss, seed_g = _grad(theta)
    g0 = float(seed_g.norm())
    gtol = float(grad_tol) if grad_tol is not None else max(grad_rtol * g0, 1e-12)
    if verbose:
        print(f"[laplace] seed: U={seed_loss:.4f}  |g|={g0:.3e}  target |g|<{gtol:.3e}")

    opt = torch.optim.Adam([theta], lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        opt, T_max=max(int(n_iters), 1), eta_min=lr / 100.0)

    best, best_theta = seed_loss, theta.detach().clone()
    bad, it, stop_reason = 0, -1, "max_iters"
    t0 = time.time()
    pbar = (tqdm(range(max(int(n_iters), 0)), desc="MAP", dynamic_ncols=True,
                 leave=True) if verbose else range(max(int(n_iters), 0)))

    for it in pbar:
        opt.zero_grad(set_to_none=False)
        lv, g = _grad(theta)
        gnorm = float(g.norm())

        if lv < best - 1e-12:
            best, best_theta = lv, theta.detach().clone()
            bad = 0
        else:
            bad += 1
        if gnorm < gtol:
            best, best_theta = lv, theta.detach().clone()
            stop_reason = "grad_converged"
            if verbose:
                pbar.set_postfix({"U": f"{lv:.4f}", "|g|": f"{gnorm:.2e}"})
            break
        if bad >= patience:
            stop_reason = f"no_improve_{patience}_iters"
            break

        theta.grad = g.detach().clone().to(dtype=theta.dtype)
        opt.step()
        sched.step()

        if verbose:
            pbar.set_postfix({"U": f"{lv:.4f}", "best": f"{best:.4f}",
                              "|g|": f"{gnorm:.2e}",
                              "lr": f"{sched.get_last_lr()[0]:.1e}",
                              "stale": f"{bad}/{patience}"})

    if verbose and hasattr(pbar, "close"): pbar.close()

    # L-BFGS
    if refine_lbfgs:
        theta_r = best_theta.clone().requires_grad_(True)
        lbfgs = torch.optim.LBFGS([theta_r], lr=1.0, max_iter=int(lbfgs_iters),
                                  history_size=20, line_search_fn="strong_wolfe",
                                  tolerance_grad=gtol, tolerance_change=1e-16)

        def closure():
            lbfgs.zero_grad(set_to_none=False)
            lv_c, g_c = _grad(theta_r)
            theta_r.grad = g_c.detach().clone().to(dtype=theta_r.dtype)
            return torch.tensor(lv_c, dtype=dtype, device=device)

        try:
            lbfgs.step(closure)
            lv_r, _ = _grad(theta_r)
            if np.isfinite(lv_r) and lv_r <= best + 1e-9:
                best, best_theta = lv_r, theta_r.detach().clone()
                stop_reason += "+lbfgs"
        except Exception as exc:
            if verbose:
                print(f"[laplace] L-BFGS skipped: {type(exc).__name__}: {str(exc)[:70]}")

    theta_map_t = best_theta.detach().clone()
    _, g_map = _grad(theta_map_t)
    gnorm_map = float(g_map.norm())
    converged = gnorm_map < gtol
    total_t = time.time() - t0

    if verbose:
        print(f"[laplace] {stop_reason} | iters={it+1} | U={best:.6f} | "
              f"|g|@MAP={gnorm_map:.3e} (target {gtol:.3e}) | {total_t:.1f}s")
    if not converged:
        print(f"[laplace] WARNING: |g| at the reported MAP is {gnorm_map:.3e}, "
              f"{gnorm_map/max(gtol,1e-30):.1f}x the target. Not a stationary "
              f"point, so the covariance below is not a posterior.")

    if hessian_mode == "autograd" and not fd_idx:
        from torch.autograd.functional import hessian as _hessian
        H = _hessian(lambda tv: _nlp_tensor(tv), theta_map_t.requires_grad_(True))
        H = 0.5 * (H + H.T)
    else:
        tnp = theta_map_t.cpu().numpy().astype(np.float64)
        H_np = np.zeros((D, D), dtype=np.float64)
        for i in range(D):
            h = _fd_step(tnp[i], fd_rel_hess)
            tp = tnp.copy(); tp[i] += h
            tm = tnp.copy(); tm[i] -= h
            _, gp = _grad(torch.tensor(tp, dtype=dtype, device=device))
            _, gm = _grad(torch.tensor(tm, dtype=dtype, device=device))
            H_np[:, i] = ((gp - gm) / (2.0 * h)).cpu().numpy()
        H = torch.tensor(0.5 * (H_np + H_np.T), dtype=dtype, device=device)

    cov, lik_eigvals, n_neg = _posterior_cov_from_hessian(H, prior_prec)

    diag = torch.clamp(torch.diag(cov), min=0.0)
    theta_map_d = {k: float(theta_map_t[i]) for i, k in enumerate(free_keys)}
    theta_std_d = {k: float(torch.sqrt(diag[i])) for i, k in enumerate(free_keys)}
    prior_std_d = {k: float(1.0 / np.sqrt(prior_prec_l[i])) for i, k in enumerate(free_keys)}

    chi2, traj_scales = None, None
    if objective == "trajectory":
        chi2, traj_scales = _trajectory_chi2(
            bayes_mm, T_list, W_list, assemble_params(theta_map_t),
            free_models, dt, D, profile_scale=traj_profile_scale)

    if verbose:
        if n_neg:
            print(f"[laplace] {n_neg}/{D} likelihood-Hessian eigenvalues negative "
                  f"(min={float(lik_eigvals.min()):.3e}); clamped to zero info, "
                  f"so those directions fall back to the prior.")
        if chi2 is not None:
            if traj_profile_scale:
                print(f"[laplace] chi2/dof = {chi2:.4f} (==1 by construction; "
                      f"fitted noise scales s^2 = "
                      f"{[f'{s:.3e}' for s in traj_scales]})")
            else:
                print(f"[laplace] chi2/dof = {chi2:.4f}  (>>1: the ODE cannot "
                      f"reproduce the coupled trajectory; <<1: W too loose, "
                      f"stds conservative)")
        print("[laplace] posterior vs prior std (ratio ~1 => the data says nothing):")
        for k in free_keys:
            print(f"[laplace]   {str(k):<24s} {theta_std_d[k]:.4f} / "
                  f"{prior_std_d[k]:.4f} = {theta_std_d[k]/prior_std_d[k]:.3f}")

    info = {
        "loss": best, "stop_reason": stop_reason, "iters": it + 1,
        "seconds": total_t, "objective": objective,
        "grad_norm_at_map": gnorm_map, "grad_norm_at_seed": g0,
        "grad_target": gtol, "converged": converged,
        "chi2_per_dof": chi2, "sensitivities": sens,
        "fd_idx": [free_keys[i] for i in fd_idx],
        "hessian": H.detach(), "hessian_eigvals": lik_eigvals.detach(),
        "n_negative_curvature": n_neg, "prior_std": prior_std_d,
        "chi2_per_dof": chi2, "traj_noise_scales": traj_scales,
        "traj_profile_scale": traj_profile_scale,
    }
    return theta_map_d, theta_std_d, cov.detach(), free_keys, info


def refit_at_map(bayes_mm, data_list, couplings, steps, params_list,
                 theta_map: Dict[Tuple[int, str], float]):
    """
    Re-run inference with the MAP parameters.
    Returns (inference_output, posterior_std_list) where posterior_std_list[m]
    has shape (T, dim).
    """
    plist = []
    for m, params in enumerate(params_list):
        d = {}
        for name, val in params.items():
            v = theta_map.get((m, name), None)
            if v is None:
                v = float(val.detach()) if torch.is_tensor(val) else float(val)
            d[name] = torch.tensor(float(v), dtype=torch.float32)
        plist.append(d)

    out = bayes_mm.inference(data_list=data_list, couplings=couplings,
                             steps=steps, params_list=plist, progress_bar=False)
    covariances_over_time = out[2]

    posterior_std_list = []
    for cov_list in covariances_over_time:
        cov_np = cov_list.detach().cpu().numpy()
        std_np = np.sqrt(np.clip(np.stack([np.diag(C) for C in cov_np], axis=0), 1e-12, None))
        posterior_std_list.append(std_np)
    return out, posterior_std_list


def delta_method(fn: Callable[[torch.Tensor], torch.Tensor],
                 theta_map: Dict[Tuple[int, str], float], cov: torch.Tensor,
                 free_keys: List[Tuple[int, str]]) -> Tuple[float, float]:
    """
    Propagate the Laplace covariance through an arbitrary scalar function of the
    free parameters, performed by autograd.
    Works for any number of free parameters and any model.
    """
    tv = torch.tensor([theta_map[k] for k in free_keys],
                      dtype=cov.dtype, requires_grad=True)
    val = fn(tv).reshape(())
    g, = torch.autograd.grad(val, tv)
    var = float(g @ cov.to(g.dtype) @ g)
    return float(val.detach()), float(np.sqrt(max(var, 0.0)))




# }}}


