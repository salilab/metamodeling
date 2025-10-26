#!/bin/bash

MATLAB_PATH="/Applications/MATLAB_R2023b.app/bin/matlab"
SCRIPT_DIR=$(dirname "$0")
METAMODEL_DIR="$SCRIPT_DIR/metamodel"

# Check if MATLAB exists
if [ ! -f "$MATLAB_PATH" ]; then
    echo "Error: MATLAB not found at $MATLAB_PATH"
    exit 1
fi

# Function to run a MATLAB script
run_matlab_script() {
    local script_name=$1
    local description=$2
    
    if [ -n "$description" ]; then
        echo "Running: $description"
    else
        echo "Running: $script_name"
    fi
    echo "Script: $script_name"
    echo "------------------------------------------"
    
    # Create a temporary MATLAB script
    temp_script=$(mktemp)
    cat > "$temp_script" <<EOF
% Add BNT to path
cd('$METAMODEL_DIR/bnt_master');
addpath(genpathKPM(pwd));

% Change to run_metamodel directory
cd('$METAMODEL_DIR/run_metamodel');

% Run the specified script
$script_name;

% Exit MATLAB
exit;
EOF
    
    # Run MATLAB in batch mode
    echo "Executing MATLAB..."
    "$MATLAB_PATH" -batch "run('$temp_script')"
    
    # Clean up
    rm -f "$temp_script"
    
    echo "Completed: $description"
    echo ""
}

# Check if a script name was provided
if [ -z "$1" ]; then
    # No argument: run all scripts
    echo "=========================================="
    echo "Running All MATLAB Metamodel Scripts"
    echo "=========================================="
    echo ""
    
    cd "$METAMODEL_DIR"
    
    # Run all MATLAB scripts for each figure
    run_matlab_script "fig3_metamodel_normal" "Generate Figure 3 data"
    run_matlab_script "fig4_metamodel_normal_test_Gpl" "Generate Figure 4 data"
    run_matlab_script "fig5_metamodel_normal_full_DGd" "Generate Figure 5 data"
    run_matlab_script "fig5_metamodel_normal_data_model" "Generate Figure 5 analysis data"
    run_matlab_script "fig6_metamodel_opt" "Generate Figure 6 optimized data"
    run_matlab_script "fig6_metamodel_non_opt" "Generate Figure 6 non-optimized data"
    run_matlab_script "fig6_metamodel_t2d_opt" "Generate Figure 6 T2D optimized data"
    run_matlab_script "fig6_metamodel_t2d_non_opt" "Generate Figure 6 T2D non-optimized data"
    
    echo "=========================================="
    echo "All MATLAB scripts completed!"
    echo "=========================================="
    echo ""
    echo "Output files are in: $METAMODEL_DIR/Output/"
    echo ""
    echo "Next step: Run Python analysis scripts"
    echo "  cd $SCRIPT_DIR/analysis"
    echo "  python Fig3.py"
    echo "  python Fig4.py"
    echo "  python Fig5.py"
    echo "  python Fig6.py"
else
    # Argument provided: run single script
    cd "$METAMODEL_DIR"
    
    run_matlab_script "$1" ""
fi

