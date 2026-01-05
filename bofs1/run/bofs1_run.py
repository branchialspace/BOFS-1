# BOFS-1 Workflow
# Usage:
# ./bofs1_run <workflow> <path/structure.cif>
# ./bofs1_run <workflow> <optimade_database> <id>

import subprocess
import sys
from pathlib import Path
import bofs1


def bofs1_test(*args):
    """BOFS-1 workflow"""
    structure_path = args[0] if len(args) == 1 else bofs1.get_structure(*args)
    name = Path(structure_path).stem
    # Serialize structure
    structure_path = bofs1.serialize_structure(structure_path)
    # QE rVV10 vc-relax structure
    relax_config = bofs1.pwx_relax_config
    relax_config['system']['input_dft'] = 'rvv10',
    relax_config['system']['vdw_corr'] = '',
    structure_path = bofs1.relax_structure(structure_path, relax_config)
    name = Path(structure_path).stem
    # spglib detect structure space group
    bofs1.spglib_structure(structure_path)
    # QE relax structure
    relax_config = bofs1.pwx_relax_config
    structure_path = bofs1.relax_structure(structure_path, relax_config)
    name = Path(structure_path).stem
    # spglib detect structure space group
    bofs1.spglib_structure(structure_path)
    # QE SCF
    scf_config = bofs1.pwx_scf_config
    np = scf_config['command'][scf_config['command'].index('-np') + 1]
    bofs1.pwx(structure_path, scf_config)
    # QE NSCF IBZ W2R
    ibz_nscf_config = bofs1.pwx_nscf_config
    ibz_nscf_config['kpts_method'] = '' # Use automatic/IBZ
    ibz_nscf_config['nosym'] = False
    ibz_nscf_config['noinv'] = False
    bofs1.pwx(structure_path, ibz_nscf_config)
    # wan2respack preprocess
    ibz_pwi = f'{name}_nscf_ibz.pwi'
    pwo = f'{name}_nscf_ibz.pwo'
    pwi = f'{name}_nscf_ibz.pwi'
    w90_config = './bofs1/wannier90/w90_configs/mlwf_config.py'
    bofs1.w90_win(pwo, pwi, w90_config, nokpts=True) # No k-points in .win for w2r .win reference
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh wan2respack_pre ./{name}/{name}.save {name} {ibz_pwi} {name}.win ./wan2respack_work', shell=True, check=True)
    subprocess.run(f'cp {name}.win {name}_ibz.win', shell=True, check=True)
    # QE NSCF W2R
    subprocess.run(f'mpirun -n {np} ./qe-7.5/bin/pw.x -nk 1 < {name}_nscf_w2r.in > {name}_nscf_w2r.out', shell=True, check=True)
    # Wannier90 preprocess W2R
    subprocess.run(f'wannier90.x -pp {name}_w2r', shell=True, check=True)
    # pw2wannier90 W2R
    pw2w90x_config = bofs1.pw2w90x_config
    pw2w90x_config['inputpp']['seedname'] = f'{name}_w2r'
    bofs1.pw2w90x(structure_path, pw2w90x_config)
    # Wannier90 W2R
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_run {name}_w2r {np}', shell=True, check=True)
    # wan2respack
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh wan2respack_post ./wan2respack_work ./respack_calc', shell=True, check=True)
    # Respack
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh respack_run ./respack_calc 1 1 {np}', shell=True, check=True)
    # QE SCF+U+V
    scf_config = bofs1.pwx_scf_config
    bofs1.pwx(structure_path, scf_config)
    # QE NSCF+U+V
    nscf_config = bofs1.pwx_nscf_config
    bofs1.pwx(structure_path, nscf_config)
    # Wannier90 preprocess
    pwo = f'{name}_nscf.pwo'
    pwi = f'{name}_nscf.pwi'
    w90_config = './bofs1/wannier90/w90_configs/mlwf_config.py'
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_preprocess {pwo} {pwi} {w90_config}', shell=True, check=True)
    # pw2wannier90
    pw2w90x_config = bofs1.pw2w90x_config
    pw2w90x_config['inputpp']['seedname'] = f'{name}'
    bofs1.pw2w90x(structure_path, pw2w90x_config)
    # Wannier90
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_run {name} {np}', shell=True, check=True)
    # QE Phonons
    phonons_config = bofs1.phx_config
    bofs1.phx(structure_path, phonons_config)
    

if __name__ == '__main__':
    workflows = {name: func for name, func in globals().items() if callable(func) and not name.startswith('_')}
    workflows[sys.argv[1]](*sys.argv[2:])















