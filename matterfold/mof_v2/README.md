Requirements (Colab Runtime):
```bash
!pip install torch_geometric
!pip install wandb
!pip install pymatgen
!pip install ase
!pip install rdkit
!gdown 1N7l6xFk58BGrvY8zZiGa6Vk_GZzsqq2U # bi_mof_129 - Bi MOFs from ARCMOF, CSDMOF, QMOF ## zip and change path, unzip here
!unzip bimofs1.zip


# ORCA
!gdown ***drive path to ORCA 6.0.1***
!chmod +x orca_6_0_1_linux_x86-64_shared_openmpi416.run
!./orca_6_0_1_linux_x86-64_shared_openmpi416.run
import os
os.environ['PATH'] = "/root/orca_6_0_1:" + os.environ['PATH']
os.environ['LD_LIBRARY_PATH'] = "/root/orca_6_0_1:" + os.environ.get('LD_LIBRARY_PATH', '')
!apt-get update && apt-get install -y tcsh
!gdown ***drive path to NBO 7***
!tar -xzvf /content/nbo7.0-bin-linux-x64.tar.gz
#!export GENEXE=/content/nbo7/bin/gennbo.i4.exe
#!export NBOEXE=/content/nbo7/bin/nbo7.i4.exe
#!export NBOBIN=/content/nbo7/bin
os.environ['PATH'] = f"/content/nbo7/bin:{os.environ['PATH']}"
os.environ['GENEXE'] = '/content/nbo7/bin/gennbo.i8.exe'
os.environ['NBOEXE'] = '/content/nbo7/bin/nbo7.i8.exe'
os.environ['NBOBIN'] = '/content/nbo7/bin'
!chmod -R +x /content/nbo7/


# QuantumESPRESSO
!wget https://www.quantum-espresso.org/rdm-download/488/v7-3-1/20c05f18ff3d167351b21c7b9043ac90/qe-7.3.1-ReleasePack.tar.gz
!tar -xf qe-7.3.1-ReleasePack.tar.gz
!apt-get install -y liblapack-dev libblas-dev libopenmpi-dev libscalapack-openmpi-dev libfftw3-dev
!mkdir -p build
!cd build
!cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DQE_FFTW_VENDOR=Internal /content/qe-7.3.1
!make -j4
# Dalcorso PAW fully-relativistic pseudopotentials
!gdown 12BcBoX8R8MSf8u0UO40UzdAfUJKpl8Qa
!unzip /content/rel_pbe.zip -d /content/rel_pbe
# ONCV fully-relativistic pseudopotentials
!git clone https://github.com/pipidog/ONCVPSP.git
!git clone https://github.com/MarioAndWario/ONCVPseudoPack.git


# GPAW
# !apt install python3-mpi4py cython3 libxc-dev gpaw-data
# !pip -q install gpaw


# MOFid failed setup
# !git clone https://github.com/snurr-group/mofid.git
# !cd /content/mofid && make init && python set_paths.py && pip install --user .
```
