#!/bin/bash
# Wrappers for RESPACK (Python2.7 env) and wan2respack (BOFS-1 Python3 env)

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../../.." && pwd )"
source "$ROOT_DIR/miniforge3/etc/profile.d/conda.sh"

# Run RESPACK in Python 2.7 env
respack_run () {
    conda activate "$ROOT_DIR/bofs1_env_py27"
    export OMP_STACKSIZE=16
    export OMP_NUM_THREADS=16
    # Forward to RESPACK binaries
    exec "$@"
}

# Run wan2respack in main BOFS-1 env
wan2respack_run () {
    conda activate "$ROOT_DIR/bofs1_env"
    exec "$ROOT_DIR/wan2respack_install/bin/wan2respack.py" "$@"
}
