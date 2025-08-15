```bash
# BOFS-1
# Ubuntu 22.04 LTS arm64

# Install
git clone https://github.com/branchialspace/BOFS-1.git
cd BOFS-1
# Set installation environment variable values in .env file
./bofs1_env.sh

# Use
cd BOFS-1
./qe_run pwx pwx_scf mofs/SIWZOO_full_n2.cif  # ./qe_run <module> <config_name> <mof_file>
