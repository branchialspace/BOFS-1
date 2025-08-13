```bash
# BOFS-1

# Installation
git clone https://github.com/branchialspace/BOFS-1.git
cd BOFS-1
# Set installation environment variable values in .env file
./bofs1_env.sh

# Use
cd BOFS-1
./qe_run <module> <config_name> <mof_file>

./qe_run pwx pwx_scf mofs/SIWZOO_full_n2.cif
