#!/usr/bin/env python3
# BOFS-1 Workflow
# Usage: ./bofs1_run.py <workflow> <structure.cif>

import subprocess
import sys
from pathlib import Path


def run(cmd):
    print(f"\n{'='*60}\n{cmd}\n{'='*60}")
    result = subprocess.run(cmd, shell=True)

def serialize(structure_path):
    path = Path(structure_path)
    if re.match(r'^\d{12}-', path.name): # check if already serialized
        return str(path)
    serial_name = f"{datetime.now():%Y%m%d%H%M}-{path.name}"
    serial_path = Path.cwd() / serial_name
    shutil.copy(path, serial_path)
    return str(serial_path)

def bofs1(structure_path):
    name = Path(structure_path).stem
    # Add Relax Step
    # QE SCF
    run(f'./qe_run pwx pwx_scf_config {structure_path}')
    # QE NSCF
    run(f'./qe_run pwx pwx_nscf_config {structure_path}')
    # Wannier90 preprocess
    pwo = f'{name}_nscf.pwo'
    pwi = f'{name}_nscf.pwi'
    config = './bofs1/wannier90/w90_configs/mlwf_config.py'
    run(f'bash ./bofs1/wannier90/w90_run.sh w90_preprocess {pwo} {pwi} {config}')
    # pw2wannier90
    run(f'./qe_run pw2w90x pw2w90x_config {structure_path}')
    # Wannier90 run
    run(f'bash ./bofs1/wannier90/w90_run.sh w90_run {name} 1')
    # wan2respack
    run(f'bash ./bofs1/respack/respack_run.sh wan2respack_pre {name} {name} {pwi} {name}.win ./wan2respack_work')
    run(f'bash ./bofs1/respack/respack_run.sh wan2respack_post ./wan2respack_work ./respack_calc')
    # RESPACK
    run(f'bash ./bofs1/respack/respack_run.sh respack_run ./respack_calc {name} input.in 16 16 1')


if __name__ == '__main__':
    workflows = {name: func for name, func in globals().items() 
                 if callable(func) and not name.startswith('_') and name not in ('run', 'serialize')}    
    structure = serialize(sys.argv[2])
    workflows[sys.argv[1]](structure)


