# Optimization of Metamodel

## Overview

This directory contains the main scripts and configuration files for running the Bayesian metamodel. The metamodel integrates three input models to simulate the glucose-stimulated insulin secretion (GSIS) system:

1. **Postprandial Model (PR)** - Postprandial glucose-insulin dynamics model
2. **Exocytosis Model (VE)** - Insulin vesicle exocytosis model  
3. **Pancreas Model (Pa)** - Pancreas-scale insulin secretion model

The metamodel is implemented using the Dynamic Bayesian Network (DBN) framework based on the MATLAB BNT toolbox.

## File Structure

### MATLAB Script Files

#### Figure Generation Scripts

- **fig3_metamodel_normal.m** - Generate Fig. 3 normal state full-time metamodel simulation
- **fig4_metamodel_normal_test_Gpl.m** - Generate Fig. 4: Selection of plasma glucose (Gpl) coupling weights by scanning different means and covariances
- **fig5_metamodel_normal_data_model.m** - Generate Fig. 5: Data-model consistency analysis with varying input model parameters
- **fig5_metamodel_normal_full_DGd.m** - Generate data for Fig. 5: Full glucose intake data simulation with different input scales
- **fig6_metamodel_opt.m** - Generate Fig. 6: Optimized normal state metamodel with updated coupling weights
- **fig6_metamodel_non_opt.m** - Generate Fig. 6: Non-optimized normal state metamodel for comparison
- **fig6_metamodel_t2d_opt.m** - Generate Fig. 6: Optimized type 2 diabetes (T2D) state metamodel
- **fig6_metamodel_t2d_non_opt.m** - Generate Fig. 6: Non-optimized T2D state metamodel for comparison

#### Testing and Validation Scripts

- **metamodel_normal_test_GVE.m** - Test glucose-vesicle exocytosis coupling variables
- **metamodel_normal_test_Gcellc.m** - Test intracellular glucose concentration
- **metamodel_normal_conflict_obs.m** - Detect conflicting observations among input models

#### Scanning and Sensitivity Analysis Scripts

- **metamodel_normal_GVE_scan_mean_as_evidence.m** - Scan GVE coupling mean as evidence
- **metamodel_normal_SVE_Gb_scan_mean_as_evidence.m** - Scan SVE and Gb means
- **metamodel_normal_SVE_Gb_scan_mean_as_evidence_full_time.m** - Full-time SVE/Gb scan

### JSON Configuration Files

Each JSON file contains model input parameters, observation data paths, and inference settings.

#### Figure Configuration Files

- **meta_normal_fig3.json** - Metamodel configuration for Fig. 3 (normal state)
- **meta_normal_fig3_t2d.json** - Metamodel configuration for Fig. 3 (T2D state)
- **meta_normal_fig4.json** - Configuration file for Fig. 4
- **meta_normal_fig5.json** - Configuration file for Fig. 5
- **meta_normal_fig6_opt.json** - Optimized version configuration for Fig. 6
- **meta_normal_fig6_non_opt.json** - Non-optimized version configuration for Fig. 6
- **meta_t2d_fig6_opt.json** - T2D optimized version configuration for Fig. 6
- **meta_t2d_fig6_non_opt.json** - T2D non-optimized version configuration for Fig. 6

#### Other Configuration Files

- **meta_normal.json** - Base normal state metamodel configuration
- **meta_normal_DGd.json** - Configuration with glucose intake data
- **meta_normal_DGd_420.json** - 420-minute timespan configuration
- **non_opt_meta_normal.json** - Simplified non-optimized version configuration
- **meta_normal copy.json** - Configuration file backup

## JSON Configuration File Format

Configuration files contain the following main sections:

```json
{
    "name": "meta_normal",
    "DataInput": {
        "// Coupling weights
        "S_ve_weight": 0.4,          // Vesicle-exocytosis coupling weight
        "G_pr_weight": 0.2,          // Plasma glucose-postprandial coupling weight  
        "G_c_weight": 0.1,           // Cellular glucose coupling weight
        "S_pa_weight": 0.2,          // Pancreas secretion coupling weight
        
        "// Observation data means and covariances
        "S_pa_obs_mean": 34.0,
        "S_cell_obs_mean": 9.32E-09,
        "G_cell_obs_mean": 2.55,
        "G_pl_obs_mean": 5.1,
        
        "// Simulation parameters
        "T": 420                      // Number of time slices (minutes)
    },
    "EvidenceInput": {
        "// Evidence data paths
        "Evidence_Scell": "../data/Scell_obs_DGd.dat",
        "Evidence_Spa": "../data/Spa_obs_DGd_avg.dat",
        "Evidence_Gcell": "../data/Gcell_obs_DGd.dat",
        "Evidence_Gpl": "../data/Gpl_obs_DGd.dat"
    }
}
```

## Key Variable Descriptions

### Coupling Variables

Coupling variables connect different input models in the metamodel framework:

- **Scell.C** - Insulin secretion rate of a single beta cell (pM min⁻¹)
- **Spa.C** - Insulin secretion rate of the pancreas (pM min⁻¹)
- **Gcell.C** - Intracellular glucose concentration (mM)
- **G.C** - Plasma glucose concentration (mM)

### Coupling PGM Graphs

The coupling probabilistic graphical model (PGM) graphs define the conditional probability distributions (CPDs) that connect surrogate model variables with coupling variables and observations:

- **P(G_pl^C | G^PR, G_pl^obs)** - Couples postprandial plasma glucose with observations
- **P(G_cell^C | G_pl^C, G_cell^obs)** - Couples plasma and intracellular glucose
- **P(S_cell^C | S^VE, S_cell^obs)** - Couples vesicle exocytosis secretion with observations
- **P(S_pa^C | S^Pa, S_pa^obs)** - Couples pancreas secretion with observations

### Weight Parameters

Weight parameters control the balance between model predictions and observational data:
- Value range: [0, 1]
- Higher weight indicates greater trust in model predictions
- Lower weight indicates greater reliance on observational data
- Weights should sum to 1 (e.g., `S_ve_weight + S_cell_obs_weight = 1`)

## Selection of Coupling Weights

### Optimization Without Reference Data

To select optimal coupling weights, the metamodel evaluates different combinations of CPD means and variances for coupling PGM graphs. The selection process balances model uncertainty reduction and model consistency.

#### Running Coupling Weight Selection

1. **For Plasma Glucose Coupling (G_pl^C)**:
   ```matlab
   fig4_metamodel_normal_test_Gpl
   ```
   This script scans through different means and variances of the CPD P(G_pl^C | G^PR, G_pl^obs).

2. **For Other Coupling Variables**:
   - Use scanning scripts to evaluate coupling parameter ranges
   - Modify JSON configuration files with candidate weight values
   - Compare model uncertainty and consistency metrics

#### Parameter Ranges for Coupling CPDs

**Mean Range**:
- Typically scanned from observed minimum to maximum values
- For glucose variables: 3.0 - 12.0 mM
- For secretion variables: 0 - 100 pM min⁻¹

**Variance Range**:
- Logarithmic scale: 10^(-4) to 10^2
- Scanned to evaluate impact of observation uncertainty
- Lower variance: stronger constraint from observations
- Higher variance: weaker constraint, more model-driven

**Weight Range**:
- Model weight: 0.0 - 1.0
- Observation weight: 1.0 - model weight
- Common tested values: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

#### Selection Criteria

The optimal coupling weights are selected by:

1. **Minimizing Model Uncertainty**: 
   - Compute Σh(X_s,upd) across all updated surrogate models
   - Lower uncertainty indicates tighter distributions

2. **Maximizing Model Consistency**: 
   - Compute η̄(X_s,upd, X_s) between surrogate and updated surrogate models
   - Higher consistency indicates better agreement

3. **Balancing Trade-off**:
   - Plot model consistency vs. model uncertainty
   - Select configuration at optimal point (e.g., "red star" in Fig. 4)

## Output Files

Results are saved in the `../Output/` directory:

- **metamodel_normal_opt.csv** - Optimized normal state metamodel results
- **metamodel_normal_non_opt.csv** - Non-optimized normal state results
- **metamodel_t2d_opt.csv** / **metamodel_t2d_non_opt.csv** - T2D state results
- **fig3_full_time_CPDs/** - Full conditional probability distributions for Fig. 3
- **fig4_test_Gpl/** - Sensitivity analysis results for Fig. 4
- **fig5_test_Gpl/** - Test results for Fig. 5

Output CSV files contain the following for all model variables at each time slice:
- Posterior mean (mu)
- Posterior covariance (Sigma)
- Posterior standard deviation (√Sigma)

## Workflow Summary

### Order of Execution

1. **Generate Surrogate Models** (Optional - outputs already provided):
   ```matlab
   cd ../surrogatemodel
   postprandial_normal
   exocytosis
   pancreas
   ```

2. **Generate Metamodel Outputs**:
   ```matlab
   cd ../run_metamodel
   % Add BNT to path first
   cd ../bnt_master
   addpath(genpathKPM(pwd))
   cd ../run_metamodel
   
   % Run MATLAB scripts for each figure
   fig3_metamodel_normal
   fig4_metamodel_normal_test_Gpl
   fig5_metamodel_normal_full_DGd
   fig6_metamodel_opt
   fig6_metamodel_non_opt
   fig6_metamodel_t2d_opt
   fig6_metamodel_t2d_non_opt
   ```

3. **Generate Figures with Python**:
   ```bash
   cd ../../analysis
   python Fig2.py  # Analyze input models
   python Fig3.py  # Metamodel full-time evolution
   python Fig4.py  # Coupling weight selection
   python Fig5.py  # Data-model consistency
   python Fig6.py  # Optimized vs non-optimized comparison
   ```

