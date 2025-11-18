#!/bin/bash
# Wrappers for RESPACK (Python2.7 env) and wan2respack (BOFS-1 Python3 env)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../../.." && pwd )"
source "$ROOT_DIR/miniforge3/etc/profile.d/conda.sh"

# Run RESPACK
respack_run () {
    local calc_dir="$1"
    local qe_save_dir="$2"
    local omp_stacksize="${3:-16}"
    local omp_num_threads="${4:-16}"
    local mpi_np="${5:-1}"
    mkdir -p "$calc_dir"
    conda activate "$ROOT_DIR/bofs1_env_py27"
    export OMP_STACKSIZE="$omp_stacksize"
    export OMP_NUM_THREADS="$omp_num_threads"
    # From Quantum ESPRESSO to RESPACK
    cp "$ROOT_DIR/RESPACK/util/qe2respack/qe2respack.py" "$calc_dir/"
    cd "$calc_dir" || return 1
    python ./qe2respack.py "$qe_save_dir"
    # Run RESPACK calculations
    mpirun -np "$mpi_np" ./calc_wannier < input.in > log.wannier
    mpirun -np "$mpi_np" ./calc_chiqw < input.in > log.chiqw
    mpirun -np "$mpi_np" ./calc_w3d < input.in > log.calc_w3d
    mpirun -np "$mpi_np" ./calc_j3d < input.in > log.calc_j3d
    # Transfer analysis
    cp "$ROOT_DIR/RESPACK/util/transfer_analysis/tr.py" "$calc_dir/"
    cd "$ROOT_DIR/RESPACK/src/transfer_analysis/" || return 1
    make
    cp calc_tr "$calc_dir/"
    cd "$calc_dir" || return 1
    python ./tr.py
    echo "RESPACK workflow completed successfully!"
}

# Run wan2respack in main BOFS-1 env
wan2respack_run () {
    conda activate "$ROOT_DIR/bofs1_env"
    exec "$ROOT_DIR/wan2respack_install/bin/wan2respack.py" "$@"
}
