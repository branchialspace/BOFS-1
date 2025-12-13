#!/bin/bash
# w90_run.sh
# Run processes for Wannier90 (BOFS-1 Python3 env)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"
source "$ROOT_DIR/miniforge3/etc/profile.d/conda.sh"

# Write win and run -pp
w90_preprocess () {
    local pwo_path="$1"
    local pwi_path="$2"
    local config_path="$3"
    conda activate "$ROOT_DIR/bofs1_env"
    local filename=$(basename -- "$pwi_path")
    local seedname="${filename%.*}"
    echo "Generating ${seedname}.win using config: ${config_path}..."
    # Run wan90_win.py to write .win file
    python "$SCRIPT_DIR/w90_win.py" "$pwo_path" "$pwi_path" "$config_path"
    # Run Wannier90 Pre-processing (-pp)
    echo "Running wannier90.x -pp on ${seedname}..."
    wannier90.x -pp "$seedname"
}

# Run wannier90.x
w90_run () {
    local input_arg="$1"
    local mpi_np="${2:-1}"
    conda activate "$ROOT_DIR/bofs1_env"
    local filename=$(basename -- "$input_arg")
    local seedname="${filename%.*}"
    echo "Running wannier90.x on ${seedname} (NP=${mpi_np})..."
    mpirun -np "$mpi_np" wannier90.x "$seedname"
}

"$@"
