#!/bin/bash
# Run processes for RESPACK (Python2.7 env) and wan2respack (BOFS-1 Python3 env)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/../../.." && pwd )"
source "$ROOT_DIR/miniforge3/etc/profile.d/conda.sh"

# Run wan2respack preprocessing on Wannier90 input
wan2respack_pre () {
    local qe_outdir="$1"
    local seedname="$2"
    local nscf_ref="$3"
    local win_ref="$4"
    local work_dir="${5:-./wan2respack_work}"
    conda activate "$ROOT_DIR/bofs1_env"
    mkdir -p "$work_dir"
    cd "$work_dir" || return 1
    cat > conf.toml <<EOF
[base]
QE_output_dir = "$qe_outdir"
seedname = "$seedname"

[pre.ref]
nscf = "$nscf_ref"
win = "$win_ref"

[pre.output]
nscf = "${seedname}_nscf_w2r.in"
win = "${seedname}.win"
EOF
    python "$ROOT_DIR/wan2respack/bin/wan2respack.py" -pp conf.toml
    cd - > /dev/null
}

# Run wan2respack Wannier90 postprocessing to prepare for RESPACK
wan2respack_post () {
    local work_dir="$1"
    local calc_dir="$2"
    conda activate "$ROOT_DIR/bofs1_env"
    cd "$work_dir" || return 1
    python "$ROOT_DIR/wan2respack/bin/wan2respack.py" conf.toml
    cp -r dir-wfn dir-wan "$calc_dir/"
    cd - > /dev/null
}

# Run RESPACK
respack_run () {
    local calc_dir="$1"
    local qe_bands_dir="$2"
    local respack_input="$3"
    local omp_stacksize="${4:-16}"
    local omp_num_threads="${5:-16}"
    local mpi_np="${6:-1}"
    mkdir -p "$calc_dir"
    conda activate "$ROOT_DIR/bofs1_env_py27"
    export OMP_STACKSIZE="$omp_stacksize"
    export OMP_NUM_THREADS="$omp_num_threads"
    cp "$ROOT_DIR/RESPACK/util/qe2respack/qe2respack.py" "$calc_dir/"
    cd "$calc_dir" || return 1
    python ./qe2respack.py "$qe_bands_dir"
    cp "$respack_input" ./input.in
    mpirun -np "$mpi_np" ./calc_wannier < input.in > log.wannier
    mpirun -np "$mpi_np" ./calc_chiqw < input.in > log.chiqw
    mpirun -np "$mpi_np" ./calc_w3d < input.in > log.calc_w3d
    mpirun -np "$mpi_np" ./calc_j3d < input.in > log.calc_j3d
    cp "$ROOT_DIR/RESPACK/util/transfer_analysis/tr.py" "$calc_dir/"
    cd "$ROOT_DIR/RESPACK/src/transfer_analysis/" || return 1
    make
    cp calc_tr "$calc_dir/"
    cd "$calc_dir" || return 1
    python ./tr.py
}

"$@"
