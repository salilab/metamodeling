# MATLAB Scripts Directory

## Overview

This directory contains all core MATLAB functions and utility scripts for building and running the Bayesian metamodel. These scripts are dependencies for the main scripts in `../run_metamodel/` directory.

## Important Notes

**These scripts are NOT meant to be run directly**. They are utility functions and modules used by the main scripts in the `../run_metamodel/` directory.

## File Categories

### 1. DBN (Dynamic Bayesian Network) Construction Functions

#### Metamodel DBN Functions
- **make_meta_dbn3.m** - Main function to create the full 3-model metamodel DBN
- **make_meta_dbn3_noobs.m** - Create metamodel DBN without observations
- **make_meta_dbn3_datamodel.m** - Create metamodel DBN for data-model consistency analysis
- **make_meta_dbn3_Gpl.m** - Create metamodel DBN focusing on plasma glucose (Gpl) coupling
- **make_meta_dbn3_Scell.m** - Create metamodel DBN focusing on single cell secretion (Scell) coupling
- **make_meta_dbn3_Gpl_Scell.m** - Create metamodel DBN focusing on both Gpl and Scell
- **make_meta_dbn3_Gcell.m** - Create metamodel DBN focusing on intracellular glucose (Gcell)
- **make_meta_dbn3_test_GVE.m** - Create metamodel DBN for testing GVE coupling

#### Surrogate Model DBN Functions

**Postprandial Model:**
- **make_postprandial_dbn.m** - Create postprandial model DBN with coupling
- **make_postprandial_dbn_pl.m** - Create postprandial model DBN for plasma level
- **make_postprandial_dbn_noobs.m** - Create postprandial DBN without observations
- **make_postprandial_dbn_nocouple.m** - Create postprandial DBN without coupling
- **make_postprandial_dbn_datamodel.m** - Create postprandial DBN for data-model analysis

**Exocytosis Model:**
- **make_exocytosis_dbn.m** - Create exocytosis model DBN with coupling
- **make_exocytosis_dbn_Scell.m** - Create exocytosis DBN focusing on Scell
- **make_exocytosis_dbn_noobs.m** - Create exocytosis DBN without observations
- **make_exocytosis_dbn_nocouple.m** - Create exocytosis DBN without coupling
- **make_exocytosis_dbn_nocouple_wDGd.m** - Create exocytosis DBN without coupling but with glucose intake data

**Pancreas Model:**
- **make_pancreas_dbn.m** - Create pancreas model DBN with coupling
- **make_pancreas_dbn_noobs.m** - Create pancreas DBN without observations
- **make_pancreas_dbn_nocouple.m** - Create pancreas DBN without coupling
- **make_pancreas_dbn_nocouple_wDGd.m** - Create pancreas DBN without coupling but with glucose intake data

#### Data Model DBN Functions
- **make_dbn_datamodel.m** - Create DBN for data-model consistency analysis
- **make_datamodel_dbn.m** - Alternative implementation for data-model DBN

### 2. DBN Factory and Merging Functions

- **DBNFactory.m** - Factory class for creating DBN instances
- **CPDFactory.m** - Factory class for creating Conditional Probability Distribution (CPD) nodes
- **merge_dbn_factories.m** - Merge multiple DBN factories into a single metamodel
- **add_bnet.m** - Add Bayesian network to the metamodel structure

### 3. Utility Functions

- **get_eclass_from_maps.m** - Get equivalence class mapping from node maps
- **get_reverse_nodes_map.m** - Get reverse mapping of nodes
- **get_valid_nodes_graph.m** - Get valid nodes and graph structure

### 4. Testing and Analysis Scripts

These scripts in the `scripts/` directory are copies or variations of main scripts:
- **fig3_metamodel_normal.m**
- **fig3_metamodel_normal_cbk.m**
- **fig4_metamodel_normal_test_Gpl.m**
- **fig5_metamodel_normal_data_model.m**
- **fig5_metamodel_normal_full_DGd.m**
- **fig6_metamodel_opt.m**
- **fig6_metamodel_non_opt.m**
- **fig6_metamodel_t2d_opt.m**
- **fig6_metamodel_t2d_non_opt.m**
- **metamodel_normal_conflict_obs.m**
- **metamodel_normal_GVE_scan_mean_as_evidence.m**
- **metamodel_normal_SVE_Gb_scan_mean_as_evidence.m**
- **metamodel_normal_SVE_Gb_scan_mean_as_evidence_full_time.m**
- **metamodel_normal_test_Gcellc.m**
- **metamodel_normal_test_GVE.m**

### 5. Standalone Surrogate Model Scripts

- **postprandial_normal.m** - Run postprandial model independently
- **exocytosis.m** - Run exocytosis model independently
- **pancreas.m** - Run pancreas model independently

**Note:** These three scripts can be run independently to generate surrogate model outputs, but they are duplicates of files in `../surrogatemodel/` directory. It's recommended to use the scripts in `../surrogatemodel/` or `../run_metamodel/` instead.

## How to Use These Scripts

These utility functions are automatically called by the main scripts in `../run_metamodel/`. To run the metamodel:

1. Navigate to the `../run_metamodel/` directory
2. Follow the instructions in `../run_metamodel/README.md`
3. Run the appropriate main script (e.g., `fig3_metamodel_normal.m`)

## Dependencies

- MATLAB R2019b or later
- BNT (Bayes Net Toolbox) - included in `../bnt_master/`
- Custom utility functions in this directory

## Related Directories

- `../run_metamodel/` - Main scripts to run for generating figures and results
- `../surrogatemodel/` - Standalone surrogate model scripts
- `../data/` - Input data files
- `../Output/` - Output results directory

## For More Information

See the main project README at the repository root for the complete workflow and usage instructions.

