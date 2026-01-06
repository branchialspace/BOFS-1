DFT environment for computing the electronic structure and optical response of MOFs

```bash
# BOFS-1
# Ubuntu 22.04 LTS x86 AMD

# Install
curl -L https://github.com/branchialspace/BOFS-1/archive/main.tar.gz | tar xz && mv BOFS-1-main BOFS-1
cd BOFS-1
bash bofs1_x86.sh

# Run bofs1
cd BOFS-1
./bofs1_run bofs1_test $MP_KEY mp-23152  # ./bofs1_run <workflow> <mp_api_key> <mp-id>
./bofs1_run bofs1_test mofs/SIWZOO_full_n2.cif  # ./bofs1_run <workflow> <path/structure.cif>

# Run QE modules
cd BOFS-1
./qe_run pwx pwx_scf_config mofs/SIWZOO_full_n2.cif  # ./qe_run <module> <config_name> <mof_file>

# Run Wannier90
cd BOFS-1
bash ./bofs1/wannier90/w90_run.sh w90_preprocess <pwo_path> <pwi_path> <config_path> # Generates win file and runs Wannier90 -pp
bash ./bofs1/wannier90/w90_run.sh w90_run <seedname> <mpi_np> # Runs wannier90.x

# Run wan2respack + RESPACK
cd BOFS-1
bash ./bofs1/respack/respack_run.sh wan2respack_pre <qe_outdir> <seedname> <nscf_ref> <win_ref> <work_dir> # For Wannier90 input
bash ./bofs1/respack/respack_run.sh wan2respack_post <work_dir> <calc_dir> # Translate wannier90 output for RESPACK
bash ./bofs1/respack/respack_run.sh respack_run <calc_dir> <qe_bands_dir> <respack_input> <omp_stacksize> <omp_num_threads> <mpi_np>

# Run PSLibrary
cd BOFS-1
bash ./bofs1/qe/pslibrary_run.sh
