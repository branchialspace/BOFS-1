# BOFS-1 Workflow
# Usage:
# ./bofs1_run <workflow> <path/structure.cif>
# ./bofs1_run <workflow> <mp_api_key> <mp-id>

import subprocess
import sys
import copy
from pathlib import Path
import bofs1


def bofs1_test(*args):
    """BOFS-1 workflow"""
    # Initialize structure
    raw_path = args[0] if len(args) == 1 else bofs1.get_structure(*args)
    serial_path = bofs1.serialize_structure(raw_path)
    sym_path = bofs1.symmetrize_structure(serial_path)
    # QE relax structure
    relax_config = copy.deepcopy(bofs1.pwx_relax_config)
    relaxed_path = bofs1.relax_structure(sym_path, relax_config)
    structure_path = bofs1.symmetrize_structure(relaxed_path)
    name = Path(structure_path).stem
    bofs1.compare_structure([sym_path, relaxed_path, structure_path])
    # QE SCF
    scf_config = copy.deepcopy(bofs1.pwx_scf_config)
    np = scf_config['command'][scf_config['command'].index('-np') + 1]
    bofs1.pwx(structure_path, scf_config)
    # QE NSCF IBZ W2R
    ibz_nscf_config = copy.deepcopy(bofs1.pwx_nscf_config)
    ibz_nscf_config['kpts_method'] = '' # Use automatic/IBZ
    ibz_nscf_config['system']['nosym'] = False
    ibz_nscf_config['system']['noinv'] = False
    bofs1.pwx(structure_path, ibz_nscf_config)
    # wan2respack preprocess
    pwo = f'{name}_nscf_ibz.pwo'
    pwi = f'{name}_nscf_ibz.pwi'
    w90_config = copy.deepcopy(bofs1.mlwf_config)
    bofs1.w90_win(pwo, pwi, w90_config, nokpts=True) # No k-points in .win for w2r .win reference
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh wan2respack_pre ./{name}/{name}.save {name} {pwi} {name}_nscf_ibz.win ./wan2respack_work', shell=True, check=True)
    # QE NSCF W2R
    subprocess.run(f'mpirun --allow-run-as-root --use-hwthread-cpus -np {np} ./qe-7.5/bin/pw.x -nk 1 < {name}_nscf_w2r.in > {name}_nscf_w2r.out', shell=True, check=True)
    # Wannier90 preprocess W2R
    subprocess.run(f'wannier90.x -pp {name}_w2r', shell=True, check=True)
    # pw2wannier90 W2R
    pw2w90x_config = copy.deepcopy(bofs1.pw2w90x_config)
    pw2w90x_config['inputpp']['seedname'] = f'{name}_w2r'
    bofs1.pw2w90x(structure_path, pw2w90x_config)
    # Wannier90 W2R
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_run {name}_w2r {np}', shell=True, check=True)
    # wan2respack
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh wan2respack_post ./wan2respack_work ./respack_calc', shell=True, check=True)
    # Respack
    subprocess.run(f'bash ./bofs1/respack/respack_run.sh respack_run ./respack_calc 1 1 {np}', shell=True, check=True)
    # QE SCF+U
    scf_config = copy.deepcopy(bofs1.pwx_scf_config)
    bofs1.pwx(structure_path, scf_config)
    # QE Phonons
    phonons_config = copy.deepcopy(bofs1.phx_config)
    bofs1.phx(structure_path, phonons_config)
    # QE NSCF+U
    nscf_config = copy.deepcopy(bofs1.pwx_nscf_config)
    nscf_config['system']['force_symmorphic'] = True
    nscf_config['nbnd_scalar'] = 4
    bofs1.pwx(structure_path, nscf_config)
    # Wannier90 preprocess
    pwo = f'{name}_nscf.pwo'
    pwi = f'{name}_nscf.pwi'
    w90_config = copy.deepcopy(bofs1.mlwf_config)
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_preprocess {pwo} {pwi} {w90_config}', shell=True, check=True)
    # pw2wannier90
    pw2w90x_config = copy.deepcopy(bofs1.pw2w90x_config)
    pw2w90x_config['inputpp']['seedname'] = f'{name}'
    bofs1.pw2w90x(structure_path, pw2w90x_config)
    # Wannier90
    subprocess.run(f'bash ./bofs1/wannier90/w90_run.sh w90_run {name} {np}', shell=True, check=True)
    # yambo

if __name__ == '__main__':
    workflows = {name: func for name, func in globals().items() if callable(func) and not name.startswith('_')}
    workflows[sys.argv[1]](*sys.argv[2:])














