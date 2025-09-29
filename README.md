DFT environment for computing the optical response of Bi-MOFs

```bash
# BOFS-1
# Ubuntu 22.04 LTS arm64

# Install
curl -L https://github.com/branchialspace/BOFS-1/archive/main.tar.gz | tar xz && mv BOFS-1-main BOFS-1
cd BOFS-1
export QE_URL=<quantum_espresso_url>
export MOF_DB_GDOWN=<mof_database_gdown_id>
export DALCORSO_GDOWN=<dalcorso_pseudopotentials_gdown_id>
bash bofs1_env.sh

# Use
cd BOFS-1
./qe_run pwx pwx_scf_config mofs/SIWZOO_full_n2.cif  # ./qe_run <module> <config_name> <mof_file>
