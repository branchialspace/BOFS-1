#!/bin/bash
# BOFS1 build environment
set -e # Exit on error
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" # ensure installation to BOFS-1 root directory
source <(sed -E 's/^([A-Za-z_][A-Za-z0-9_]*)=(.*)$/export \1=\${\1:-\2}/' .env) # export installation .env variables. parameter expansion defaults to precedent
# miniforge3
curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
bash Miniforge3.sh -b -p ./miniforge3
rm -f Miniforge3.sh
source ./miniforge3/etc/profile.d/conda.sh
source ./miniforge3/etc/profile.d/mamba.sh
# Create mamba venv
mamba create -y -p ./bofs1_env python=3.10
conda activate ./bofs1_env
# dependencies
mamba install -y -c conda-forge cmake git wget unzip openmpi openmpi-mpicc fftw lapack blas scalapack
pip install numpy==1.26.4
pip install torch_geometric
pip install wandb
pip install pymatgen
pip install ase
pip install rdkit
pip install mendeleev
pip install gdown
gdown $MOF_DB_GDOWN # 369 Bi MOFs from ARCMOF, CSDMOF, QMOF, MOSAEC-DB
mkdir -p mofs
unzip bimofs2.zip -d mofs
# QuantumESPRESSO
wget $QE_URL
tar -xzf qe-7.4.1-ReleasePack.tar.gz
cmake -DCMAKE_C_COMPILER=$CONDA_PREFIX/bin/mpicc -DCMAKE_Fortran_COMPILER=$CONDA_PREFIX/bin/mpif90 -DQE_FFTW_VENDOR=Internal -DQE_ENABLE_OPENMP=ON qe-7.4.1
make -j4
make ld1
mkdir -p qe-7.4.1/bin
cp bin/ld1.x qe-7.4.1/bin/
# Dalcorso fully-relativistic pseudopotentials
gdown $DALCORSO_GDOWN
unzip dalcorso_rel_pbe.zip
# ONCV fully-relativistic pseudopotentials repositories
git clone https://github.com/pipidog/ONCVPSP.git
git clone https://github.com/MarioAndWario/ONCVPseudoPack.git
# BOFS1 QE runner venv wrapper
cat > qe_run << 'EOF'
#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$SCRIPT_DIR/miniforge3/etc/profile.d/conda.sh"
conda activate "$SCRIPT_DIR/bofs1_env"
exec python "$SCRIPT_DIR/bofs1/qe/qe_run.py" "$@"
EOF
chmod +x qe_run

echo "BOFS1 environment built successfully"
