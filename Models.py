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


def ODE_eqs_model_PXR(C, t, k, **kwargs):
    
    C_PXR, C_RNA, C_RIF, C_CYP = C
    CYP_0, RNA_0, EC50_mean, EC50_stdev, k_inact, k_rnadeg,\
    k_cypdeg, p, q = k
    #EC50 = torch.normal(EC50_mean, EC50_stdev)
    EC50 = EC50_mean
    CYP_prime = C_CYP/CYP_0
    RNA_prime = C_RNA
    dC_PXR = (1+p)/(1+p*CYP_prime)*(C_RIF)/(EC50+C_RIF)\
        -k_inact*C_PXR
    dC_RNA_prime = k_rnadeg*(1+q*C_PXR-RNA_prime)
    d_CYP_prime = k_cypdeg*(RNA_prime-CYP_prime)
    dC_RNA = dC_RNA_prime
    d_CYP = d_CYP_prime*CYP_0
    dC_RIF = 0*d_CYP

    return np.array([dC_PXR, dC_RNA, dC_RIF, d_CYP],dtype=object)

def ODE_eqs_model_MM(C, t, k, **kwargs):
    r"""
    Model for computing clearance and metabolism

    """
    C_substrate, C_CYP, CL = C
    C3A4_in, C3A4_sub, CYP_0, kcat, Km, Ki = k

    isInh = 0.00001*(C3A4_in)
    isSub = (C3A4_sub)
    kcatByKm = kcat/Km
    #print(kcatByKm)
    CL = isSub*kcatByKm*(1/(1+isInh*C_substrate/Ki))*C_CYP
    #print(CL)
    dC_substrate = -CL*C_substrate
    dC_CYP = 0*C_CYP
    dCL = 0*dC_CYP
    return np.array([dC_substrate,dC_CYP,dCL],dtype=object)

def ODE_eqs_model_ADMETAI(C, t, k, **kwargs):
    r"""
    Model for computing clearance

    """
    C_substrate, CL = C
    k_CL = k[0]
    CL = k_CL
    dC_substrate = -CL*C_substrate
    dCL = 0*dC_substrate
    return np.array([dC_substrate,dCL],dtype=object)


def ODE_eqs_model_human(C, t, k, **kwargs):
    r"""
    Model A ODE System

    """
    C_p, C_lung, C_cell, C_case, CL = C
    k1, k2, k3, k4, k5, k6, k7, VC_human = k
    dC_p = -k1 * C_p - CL/VC_human*C_p
    dC_lung = k2*k3*(C_p-C_lung)
    dC_cell = k4*k5*(C_p-C_cell)
    dC_case = k6*k7*(C_p-C_case)
    dCL = 0*dC_p
    return np.array([dC_p, dC_lung, dC_cell, dC_case, dCL],dtype=object)


def ODE_eqs_model_rabbit(C, t, k, **kwargs):
    r"""
    Model B ODE System

    """
    C_p, P, C_lung, C_cell, C_case = C
    k1, k2, k3, k4, k5, k6, k7, Q, V_C, V_P = k
    dC_p = -k1 * C_p - (Q/V_C)*C_p + (Q/V_P)*P
    dP =  (Q/V_C)*C_p - (Q/V_P)*P
    dC_lung = k2*k3*(C_p-C_lung)
    dC_cell = k4*k5*(C_p-C_cell)
    dC_case = k6*k7*(C_p-C_case)
    return np.array([dC_p, dP, dC_lung, dC_cell, dC_case],dtype=object)

