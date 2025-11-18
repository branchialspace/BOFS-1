DFT environment for computing the electronic structure and optical response of MOFs

```bash
# BOFS-1
# Ubuntu 22.04 LTS x86 AMD

# Install
curl -L https://github.com/branchialspace/BOFS-1/archive/main.tar.gz | tar xz && mv BOFS-1-main BOFS-1
cd BOFS-1
bash bofs1_x86.sh

# Run QE modules
cd BOFS-1
./qe_run pwx pwx_scf_config mofs/SIWZOO_full_n2.cif  # ./qe_run <module> <config_name> <mof_file>

# Run PSLibrary
cd BOFS-1
bash ./bofs1/qe/pslibrary_run.sh

# Run RESPACK
cd BOFS-1
bash ./bofs1/respack/respack_run.sh <new_calc_dir> <qe_bands_out_dir> <omp_stacksize> <omp_num_threads> <mpi_np>
