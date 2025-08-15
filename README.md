```bash
# BOFS-1
# Ubuntu 22.04 LTS arm64

# Install
git clone https://github.com/branchialspace/BOFS-1.git
cd BOFS-1
export QE_URL=<quantum_espresso_url>
export MOF_DB_GDOWN=<mof_database_gdown_id>
export DALCORSO_GDOWN=<dalcorso_pseudopotentials_gdown_id>
./bofs1_env.sh

# Use
cd BOFS-1
./qe_run pwx pwx_scf mofs/SIWZOO_full_n2.cif  # ./qe_run <module> <config_name> <mof_file>
