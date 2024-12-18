Requirements (Colab Runtime):
```bash
!pip install rdkit
!pip install ase
!pip install pymatgen
!git clone https://github.com/branchialspace/matterfold.git

!gdown ***drive path to ORCA 6.0.1***
!chmod +x orca_6_0_1_linux_x86-64_shared_openmpi416.run
!./orca_6_0_1_linux_x86-64_shared_openmpi416.run
import os
os.environ['PATH'] = "/root/orca_6_0_1:" + os.environ['PATH']
os.environ['LD_LIBRARY_PATH'] = "/root/orca_6_0_1:" + os.environ.get('LD_LIBRARY_PATH', '')

!apt-get update && apt-get install -y tcsh
!gdown ***drive path to NBO7***
!tar -xzvf /content/nbo7.0-bin-linux-x64.tar.gz
#!export GENEXE=/content/nbo7/bin/gennbo.i4.exe
#!export NBOEXE=/content/nbo7/bin/nbo7.i4.exe
#!export NBOBIN=/content/nbo7/bin
os.environ['PATH'] = f"/content/nbo7/bin:{os.environ['PATH']}"
os.environ['GENEXE'] = '/content/nbo7/bin/gennbo.i8.exe'
os.environ['NBOEXE'] = '/content/nbo7/bin/nbo7.i8.exe'
os.environ['NBOBIN'] = '/content/nbo7/bin'
!chmod -R +x /content/nbo7/

# !apt install python3-mpi4py cython3 libxc-dev gpaw-data
# !pip -q install gpaw
```
Example Bi2C4H2O4 MOF Nanoparticle:
![](bi2_fumarate_nanoparticle.png)
