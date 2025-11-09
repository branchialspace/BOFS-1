#!/bin/bash
# BOFS1 build environment
# Ubuntu 22.04 LTS x86 AMD
set -exo pipefail # Exit on error
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR" # ensure installation to BOFS-1 root directory
# miniforge3
curl -L -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3.sh -b -p ./miniforge3
rm -f Miniforge3.sh
source ./miniforge3/etc/profile.d/conda.sh
export MAMBA_ROOT_PREFIX="$SCRIPT_DIR/miniforge3"
# Create mamba venv
mamba create -y -p ./bofs1_env python=3.10
conda activate ./bofs1_env
# dependencies
mamba install -y -c conda-forge cmake ninja git wget unzip openmpi openmpi-mpicc gfortran binutils
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
# AMD Optimizing CPU Libraries
wget https://download.amd.com/developer/eula/aocl/aocl-5-1/aocl-linux-gcc-5.1.0_1_amd64.deb
sudo dpkg -i aocl-linux-gcc-5.1.0_1_amd64.deb
export AOCL_DIR=/opt/AMD/aocl/aocl-linux-gcc-5.1.0/gcc
export AOCL_LIB=$AOCL_DIR/lib_LP64
export LD_LIBRARY_PATH=$AOCL_LIB:$LD_LIBRARY_PATH
# Quantum ESPRESSO
git clone https://gitlab.com/QEF/q-e.git qe-7.5 && (cd qe-7.5 && git checkout -b qe-7.5-pinned 770a0b2d12928a67048e2f3da8d10d057e52179e)
cmake -G Ninja \
  -S qe-7.5 \
  -B build_qe \
  -DCMAKE_INSTALL_PREFIX=$SCRIPT_DIR/qe-7.5 \
  -DCMAKE_C_COMPILER=$CONDA_PREFIX/bin/mpicc \
  -DCMAKE_Fortran_COMPILER=$CONDA_PREFIX/bin/mpif90 \
  -DCMAKE_C_FLAGS="-O3 -march=znver3 -mtune=znver3" \
  -DCMAKE_Fortran_FLAGS="-O3 -march=znver3 -mtune=znver3 -fallow-argument-mismatch" \
  -DCMAKE_EXE_LINKER_FLAGS="-no-pie" \
  -DQE_ENABLE_OPENMP=ON \
  -DQE_ENABLE_SCALAPACK=ON \
  -DQE_FFTW_VENDOR=FFTW3 \
  -DBLAS_LIBRARIES="$AOCL_LIB/libblis-mt.a;-lm;-lpthread;-lgfortran" \
  -DLAPACK_LIBRARIES="$AOCL_LIB/libflame.so;-lm;-lpthread;-lgfortran" \
  -DFFTW_LIBRARIES="$AOCL_LIB/libfftw3.a;$AOCL_LIB/libfftmpi.a" \
  -DSCALAPACK_LIBRARIES="$AOCL_LIB/libscalapack.a"
ninja -C build_qe
ninja -C build_qe install
rm -rf "$SCRIPT_DIR/build_qe"
echo "export PATH=$SCRIPT_DIR/qe-7.5/bin:\$PATH" >> "$CONDA_PREFIX/etc/conda/activate.d/qe_path.sh"
echo "export LD_LIBRARY_PATH=$SCRIPT_DIR/qe-7.5/lib:/opt/AMD/aocl/aocl-linux-gcc-5.1.0/gcc/lib_LP64:\$LD_LIBRARY_PATH" >> "$CONDA_PREFIX/etc/conda/activate.d/qe_path.sh"
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
# Project-level permissions
chmod -R u+rwX,go+rwX "$SCRIPT_DIR"

echo "built BOFS1 environment"
