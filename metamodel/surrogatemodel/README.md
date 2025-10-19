# Surrogate Models

## Overview

This directory contains the surrogate models used in the Bayesian metamodel. Each surrogate model independently describes a specific biological scale of the glucose-stimulated insulin secretion (GSIS) system and is converted to a probabilistic Dynamic Bayesian Network (DBN) representation.

## Model List

### 1. Postprandial Model (PR)
**Description**: Describes the system dynamics of postprandial glucose and insulin

**Primary Variables**:
- `G` - Plasma glucose concentration (mM)
- `DG` - Glucose concentration above basal level (mM)
- `DGd` - Rate of glucose intake from food digestion (mM min⁻¹)
- `Gb` - Basal glucose concentration (mM)
- `I` - Plasma insulin concentration (pM)
- `S` - Pancreatic insulin secretion rate (pM min⁻¹)
- `Sb` - Basal insulin secretion rate (pM min⁻¹)
- `Y` - Beta cell new insulin supply (pM min⁻¹)

**Key Parameters**:
- `alpha` - Delay between glucose signal and insulin secretion (min⁻¹)
- `beta` - Pancreatic responsiveness to glucose (pM min⁻¹ mM⁻¹)
- `gamma` - Transport rate between portal vein and liver (min⁻¹)
- `k1` - Coefficient for insulin-mediated glucose reduction (min⁻¹)
- `k2` - Glucose self-reduction coefficient (min⁻¹)
- `k3` - Elevated glucose self-reduction coefficient (min⁻¹)
- `k4` - Insulin secretion degradation coefficient
- `K` - Pancreatic responsiveness to glucose change rate (pmol L⁻¹ mM⁻¹)

### 2. Exocytosis Model (VE)
**Description**: Describes the insulin vesicle exocytosis process in beta cells

**Primary Variables**:
- `G` - Intracellular glucose concentration (mM)
- `kt` - Effective rate of vesicle trafficking towards cellular periphery (m s⁻¹)
- `Npatch` - Number of activation patches per vesicle
- `Nvesicle` - Number of insulin vesicles in one beta cell
- `Ninsulin` - Amount of insulin in one vesicle (pmol)
- `Rcell` - Beta cell radius (μm)
- `Dvesicle` - Diffusion coefficient of insulin vesicles in beta cell (Å² fs⁻¹)
- `S` - Insulin secretion rate of one beta cell (pM min⁻¹)

**Key Parameters**:
- `alpha` - Correlation between secretion rate and transport force coefficient (m s⁻¹ pM⁻¹ min)
- `beta` - Correlation between secretion rate and vesicle number (pM⁻¹ min)
- `kG` - Coefficient for glucose-stimulated insulin secretion (pM min⁻¹ mM⁻¹)
- `kp` - Coefficient for activation patches accelerating secretion (pM min⁻¹)
- `kinsulin` - Coefficient for insulin amount in vesicles determining secretion (pM min⁻¹)
- `kD` - Coefficient for vesicle diffusion promoting secretion (pM min⁻¹ Å⁻² fs)
- `kR` - Coefficient for cell radius reducing secretion (pM min⁻¹ μm⁻¹)

### 3. Pancreas Model (Pa)
**Description**: Scales single-cell insulin secretion to whole pancreas level

**Primary Variables**:
- `Scell` - Insulin secretion rate of a single beta cell (pM min⁻¹)
- `Sislet` - Insulin secretion rate of a single islet (pM min⁻¹)
- `Spancreas` - Insulin secretion rate of the pancreas (pM min⁻¹)

**Key Parameters**:
- `Nc` - Number of beta cells in an islet (~1140)
- `Ni` - Number of islets in a pancreas (~3.2×10⁶)

## File Structure

### MATLAB Scripts

- **postprandial_normal.m** - Run independent simulation of postprandial model
- **exocytosis.m** - Run independent simulation of exocytosis model
- **pancreas.m** - Run independent simulation of pancreas model

### JSON Configuration Files

#### Postprandial Model Configurations

- **postprandial_normal_s1.json** - Normal state scenario 1
- **postprandial_normal_s2.json** - Normal state scenario 2
- **postprandial_non_opt_normal.json** - Normal state non-optimized version
- **postprandial_non_opt_t2d.json** - T2D state non-optimized version
- **postprandial_t2d_s1.json** - T2D state scenario 1
- **postprandial_t2d_s2.json** - T2D state scenario 2

#### Other Model Configurations

- **exocytosis.json** - Exocytosis model configuration
- **pancreas.json** - Pancreas model configuration

## JSON Configuration File Format

Each configuration file contains three main sections:

### 1. DataInput - Model Variables and Parameters

```json
{
    "name": "model_name",
    "DataInput": {
        "_comment": "Model variables",
        "variable_mean": 5.1,
        "variable_cov": 0.01,
        
        "_comment": "Model parameters",  
        "parameter": 0.05,
        
        "_comment": "Inference parameters",
        "dt": 1,
        "cov_scale": 1E-2,
        "T": 420
    }
}
```

### 2. EvidenceInput - Evidence Data

Points to data files containing observational evidence.

### 3. EvidenceDG/EvidenceScell - Input Data

Points to data files containing model driving inputs (e.g., glucose intake rate).

## Output Files

Results are saved in the `../Output/` directory:

- **postprandial_prior_normal_wDGd.csv** - Postprandial model output
- **exocytosis_prior_wDGd.csv** - Exocytosis model output
- **pancreas_prior_wDGd.csv** - Pancreas model output

Output includes the following for each variable at all time slices:
- Posterior mean (mu)
- Posterior covariance (Sigma)
- Posterior standard deviation (sqrt(Sigma))

