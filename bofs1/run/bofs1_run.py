#!/usr/bin/env python3
# BOFS-1 Workflow
# Usage: ./bofs1_run.py <structure.cif>

import subprocess
import sys
from pathlib import Path


def run(cmd):
    print(f"\n{'='*60}\n{cmd}\n{'='*60}")
    result = subprocess.run(cmd, shell=True)

def bofs1(structure_path):
    """Full workflow: SCF -> NSCF -> Wannier90 -> RESPACK"""
    name = Path(structure_path).stem
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
    bofs1(sys.argv[1])
