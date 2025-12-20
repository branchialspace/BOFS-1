# BOFS-1 Workflow
# Usage: ./bofs1_run.py <workflow> <structure.cif>

import subprocess
import sys
import shutil
import re
from pathlib import Path
from datetime import datetime
import bofs1


def serialize(structure_path):
    """Add timestamp to structure filename"""
    path = Path(structure_path)
    if re.match(r'^\d{12}-', path.name):  # check if already serialized
        return str(path)
    serial_name = f"{datetime.now():%Y%m%d%H%M}-{path.name}"
    serial_path = Path.cwd() / serial_name
    shutil.copy(path, serial_path)
    return str(serial_path)

def bofs1_run(structure_path):
    """BOFS-1 workflow"""
    name = Path(structure_path).stem
    # QE SCF
    scf_config = bofs1.pwx_scf_config
    bofs1.pwx(structure_path, scf_config)
    # QE NSCF
    nscf_config = bofs1.pwx_nscf_config
    bofs1.pwx(structure_path, nscf_config)
    # Wannier90 preprocess
    pwo = f'{name}_nscf.pwo'
    pwi = f'{name}_nscf.pwi'
    w90_config = './bofs1/wannier90/w90_configs/mlwf_config.py'
    subprocess.run(
        f'bash ./bofs1/wannier90/w90_run.sh w90_preprocess {pwo} {pwi} {w90_config}', 
        shell=True, check=True)
    # pw2wannier90
    bofs1.pw2w90x(structure_path, bofs1.pw2w90x_config)
    subprocess.run(
        # Wannier90
        f'bash ./bofs1/wannier90/w90_run.sh w90_run {name} 1 && '
        # wan2respack
        f'bash ./bofs1/respack/respack_run.sh wan2respack_pre {name} {name} {pwi} {name}.win ./wan2respack_work && '
        f'bash ./bofs1/respack/respack_run.sh wan2respack_post ./wan2respack_work ./respack_calc && '
        # Respack
        f'bash ./bofs1/respack/respack_run.sh respack_run ./respack_calc {name} input.in 16 16 1',
        shell=True, check=True)

if __name__ == '__main__':
    workflows = {name: func for name, func in globals().items() if callable(func) and not name.startswith('_') and name != 'serialize'}
    structure = serialize(sys.argv[2])
    workflows[sys.argv[1]](structure)
