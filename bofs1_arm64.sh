#!/bin/bash
# BOFS1 build environment
# Ubuntu 22.04 LTS arm64
set -e # Exit on error
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" # ensure installation to BOFS-1 root directory
source <(sed -E 's/^([A-Za-z_][A-Za-z0-9_]*)=(.*)$/export \1=\${\1:-\2}/' .env) # export installation .env variables. parameter expansion defaults to precedent
# miniforge3
curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"
bash Miniforge3.sh -b -p ./miniforge3
rm -f Miniforge3.sh
source ./miniforge3/etc/profile.d/conda.sh
export MAMBA_ROOT_PREFIX="$SCRIPT_DIR/miniforge3"
# Create mamba venv
mamba create -y -p ./bofs1_env python=3.10
conda activate ./bofs1_env
# dependencies
mamba install -y -c conda-forge cmake ninja git wget unzip openmpi openmpi-mpicc fftw lapack blas scalapack gfortran binutils
pip install numpy==1.26.4
pip install torch_geometric
pip install wandb
pip install pymatgen
pip install ase
pip install rdkit
pip install mendeleev
pip install gdown
gdown 1p4Pjl8_nrV4lYY_vIZ6dn4tseQ7iTY1v # 369 Bi MOFs from ARCMOF, CSDMOF, QMOF, MOSAEC-DB
mkdir -p mofs
unzip bimofs2.zip -d mofs
# ARM Performance Libraries
wget https://developer.arm.com/-/cdn-downloads/permalink/Arm-Performance-Libraries/Version_25.07.1/arm-performance-libraries_25.07.1_deb_gcc.tar
tar -xvf arm-performance-libraries_25.07.1_deb_gcc.tar
sudo bash ./arm-performance-libraries_25.07.1_deb/arm-performance-libraries_25.07.1_deb.sh --accept
export ARMPL_DIR=/opt/arm/armpl_25.07.1_gcc
# Quantum ESPRESSO
git clone https://gitlab.com/QEF/q-e.git qe-7.5 && (cd qe-7.5 && git checkout -b qe-7.5-pinned 770a0b2d12928a67048e2f3da8d10d057e52179e)
cmake -G Ninja \
  -DCMAKE_C_COMPILER=$CONDA_PREFIX/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=$CONDA_PREFIX/bin/mpif90 \
  -DQE_FFTW_VENDOR=Internal \
  -DQE_ENABLE_OPENMP=ON \
  -DBLAS_LIBRARIES="$ARMPL_DIR/lib/libarmpl_lp64_mp.so;-lgfortran;-lm;-lpthread" \
  -DLAPACK_LIBRARIES="$ARMPL_DIR/lib/libarmpl_lp64_mp.so;-lgfortran;-lm;-lpthread" \
  qe-7.5
ninja
ninja ld1
mkdir -p qe-7.5/bin
cp bin/ld1.x qe-7.5/bin/
# Dalcorso fully-relativistic pseudopotentials
git clone https://github.com/dalcorso/pslibrary.git
sed -i "s|PWDIR='/path_to_quantum_espresso/'|PWDIR='../../qe-7.5'|" ./pslibrary/QE_path
bash ./bofs1/qe/pslibrary_run.sh
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

echo "built BOFS1 environment"
