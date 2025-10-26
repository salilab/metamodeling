## A Bayesian approach for optimizing metamodeling under uncertainty: Application to complex biological systems

### Summary
Bayesian metamodeling divides and conquers this task by integrating a collection of input models. It proceeds through three stages: (i) convert input models into probabilistic surrogate models; (ii) couple surrogate models; and (iii) update surrogate models and input models via backpropagation. Comprehensive understanding and quantitative assessment of uncertainty in Bayesian metamodeling are crucial for its proper interpretation, although challenging due to its inherent high-complexity. We analyze the propagation of uncertainty across metamodeling stages using both an analytical example of two Gaussian time-independent input models and a numerical example of three time-dependent input models describing glucose-stimulated insulin secretion in pancreatic $\beta$-cells. We elucidate criteria in the selection of surrogate models and coupling strategies, as well as identification and resolution of potential conflicts among the input models. The optimized metamodel demonstrates better alignment with experimental data compared to the non-optimized one. 

### Project Structure

```
metamodeling/
├── README.md                          # This file - project overview and quick start guide
├── metamodel/                         # MATLAB implementation of the metamodel
│   ├── bnt_master/                    # BNT (Bayes Net Toolbox) by Kevin Murphy
│   ├── scripts/                       # Core MATLAB utility functions (see scripts/README.md)
│   ├── run_metamodel/                 # Main MATLAB scripts to generate figures (see run_metamodel/README.md)
│   ├── surrogatemodel/                # Independent surrogate model scripts (see surrogatemodel/README.md)
│   ├── data/                          # Input data files (observations, parameters, etc.)
│   └── Output/                        # Output files from MATLAB scripts (CSV format)
└── analysis/                          # Python scripts for analysis and figure generation (see analysis/README.md)
    ├── analyze_metamodel.py           # Core analysis utility functions
    ├── statistics_basic.py            # Basic statistical functions
    ├── Fig2.py                        # Generate Figure 2
    ├── Fig3.py                        # Generate Figure 3
    ├── Fig4.py                        # Generate Figure 4
    ├── Fig5.py                        # Generate Figure 5
    ├── Fig6.py                        # Generate Figure 6
    └── input_pr/                      # Input observational data for validation
```

### Directory Overview

- **`metamodel/`** - Contains MATLAB scripts for constructing the glucose-stimulated-insulin-secretion (GSIS) metamodel with three input models, using the BNT package by Kevin Murphy [http://github.com/bayesnet/bnt](http://github.com/bayesnet/bnt).
  - **`bnt_master/`** - Bayes Net Toolbox for MATLAB
  - **`scripts/`** - Core utility functions and DBN construction functions
  - **`run_metamodel/`** - Main scripts to run for generating figure data
  - **`surrogatemodel/`** - Scripts for running individual surrogate models
  - **`data/`** - Input data files (observations, model parameters, etc.)
  - **`Output/`** - Generated output files from MATLAB scripts

- **`analysis/`** - Contains Python scripts for metamodel analysis and generating figures in the manuscript.

### Workflow Overview

The project follows a two-stage workflow:

**Stage 1: MATLAB - Generate Metamodel Outputs**
1. Run MATLAB scripts in `metamodel/run_metamodel/` to generate simulation data
2. Outputs are saved as CSV files in `metamodel/Output/`

**Stage 2: Python - Analyze and Visualize**
1. Run Python scripts in `analysis/` to process MATLAB outputs
2. Generate publication-quality figures


### Quick Start Guide

#### Prerequisites

**MATLAB Requirements:**
- MATLAB R2019b or later
- BNT toolbox (included in `metamodel/bnt_master/`)

**Python Requirements:**
- Python 3.7 or later
- Required packages: `numpy`, `scipy`, `pandas`, `matplotlib`

Install Python dependencies:
```bash
pip install numpy scipy pandas matplotlib
```

#### Step 1: Setup MATLAB Path

Open MATLAB and add the BNT toolbox to the path:
```matlab
cd /path/to/metamodeling/metamodel/bnt_master
addpath(genpathKPM(pwd))
```

#### Step 2: Run MATLAB Scripts

You can run MATLAB scripts in two ways:

**Option A: Using Shell Scripts**

Run all scripts at once:
```bash
./run_matlab.sh
```

Or run individual scripts:
```bash
./run_matlab.sh fig3_metamodel_normal
./run_matlab.sh fig4_metamodel_normal_test_Gpl
./run_matlab.sh fig5_metamodel_normal_full_DGd
./run_matlab.sh fig6_metamodel_opt
./run_matlab.sh fig6_metamodel_non_opt
./run_matlab.sh fig6_metamodel_t2d_opt
./run_matlab.sh fig6_metamodel_t2d_non_opt
```

**Option B: Manual MATLAB Execution**

Open MATLAB and run scripts manually:
```matlab
cd metamodel/bnt_master
addpath(genpathKPM(pwd))
cd ../run_metamodel

% Generate Figure 3 data
fig3_metamodel_normal

% Generate Figure 4 data
fig4_metamodel_normal_test_Gpl

% Generate Figure 5 data
fig5_metamodel_normal_full_DGd

% Generate Figure 6 data
fig6_metamodel_opt
fig6_metamodel_non_opt
fig6_metamodel_t2d_opt
fig6_metamodel_t2d_non_opt
```

#### Step 3: Run Python Scripts

```bash
cd /path/to/metamodeling/analysis

# Generate Figure 2
python Fig2.py

# Generate Figure 3
python Fig3.py

# Generate Figure 4
python Fig4.py

# Generate Figure 5
python Fig5.py

# Generate Figure 6
python Fig6.py
```
