#!/bin/bash
# BOFS1 build environment
set -e # Exit on error
python3 -m venv bofs1_env
source bofs1_env/bin/activate
sudo apt-get update
sudo apt-get install -y python3-pip
pip install numpy==1.26.4
pip install torch_geometric
pip install wandb
pip install pymatgen
pip install ase
pip install rdkit
pip install mendeleev
pip install gdown
gdown ... # 369 Bi MOFs from ARCMOF, CSDMOF, QMOF, MOSAEC-DB
mkdir -p mofs
unzip bimofs2.zip -d mofs
# QuantumESPRESSO
wget ...
tar -xzf qe-7.4.1-ReleasePack.tar.gz
sudo apt-get update
sudo apt-get install -y liblapack-dev libblas-dev libopenmpi-dev libscalapack-openmpi-dev libfftw3-dev
cmake -DCMAKE_C_COMPILER=mpicc -DQE_FFTW_VENDOR=Internal -DQE_ENABLE_OPENMP=ON qe-7.4.1
make -j4
make ld1
mkdir -p qe-7.4.1/bin
cp bin/ld1.x qe-7.4.1/bin/
# Dalcorso fully-relativistic pseudopotentials
gdown ...
unzip dalcorso_paw_fr.zip -d ..
# ONCV fully-relativistic pseudopotentials repositories
git clone https://github.com/pipidog/ONCVPSP.git
git clone https://github.com/MarioAndWario/ONCVPseudoPack.git

echo "Environment setup completed successfully!"
