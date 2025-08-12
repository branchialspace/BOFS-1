```bash
# BOFS-1 Installation

# Clone
git clone https://github.com/branchialspace/BOFS-1.git

# Set installation .env variables
export QE_URL=<quantum_espresso_url>
export MOF_DB_GDOWN=<mof_database_gdown_id>
export DALCORSO_GDOWN=<dalcorso_pseudopotentials_gdown_id>

# Run installation
./BOFS-1/bofs1_env.sh
