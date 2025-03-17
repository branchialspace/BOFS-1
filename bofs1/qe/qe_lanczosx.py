# Run QuantumESPRESSO turbo_lanczos.x for TDDFT-Lanczos calculations

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
from ase.io import read


def qe_lanczosx(
    structure_path,
    config
):
    """
    Run QE turbo_lanczos.x calculation for TDDFT-Lanczos response properties.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_lr_input(config, input_filename):
        """
        Write the QE turbo_lanczos.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # LR_INPUT namelist (always required)
            f.write('&LR_INPUT\n')
            for key, value in config['lr_input'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')    
            # LR_CONTROL namelist
            f.write('\n&LR_CONTROL\n')
            for key, value in config['lr_control'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')
            # LR_POST namelist (optional, only if charge_response is set to 1)
            if config['lr_control'].get('charge_response', 0) == 1:
                f.write('\n&LR_POST\n')
                for key, value in config['lr_post'].items():
                    if isinstance(value, bool):
                        val = '.true.' if value else '.false.'
                    elif isinstance(value, str):
                        val = f"'{value}'"
                    else:
                        val = value
                    f.write(f"  {key} = {val}\n")
                f.write('/\n')
    
    # Args
    structure = read(structure_path)  # ASE Atoms object
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "turbo_lanczos"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['lr_input']['prefix'] = structure_name
    config['lr_input']['outdir'] = structure_name
    os.makedirs(structure_name, exist_ok=True)
    # Write QE turbo_lanczos.x input file
    write_lr_input(config, f"{run_name}.tdfi")
    # Subprocess run
    try:
        with open(f"{run_name}.tdfo", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.tdfi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("TDDFT-Lanczos calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running TDDFT-Lanczos calculation: {cpe}")
        try:
            with open(f"{run_name}.tdfo", 'r') as f_out:
                print("\nTDDFT-Lanczos Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


lanczos_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_lanczos.x'],
    'lr_input': {
        'wfcdir': './',
        'restart_step': 500,     # Write restart files every restart_step iterations
        'lr_verbosity': 1
    },
    'lr_control': {
        'itermax': 500,          # Number of Lanczos iterations
        'ipol': 4,               # 1=X, 2=Y, 3=Z, 4=all polarizations
        'n_ipol': 3,             # Number of zeta coefficients to calculate
        'ltammd': False,         # Tamm-Dancoff approximation
        'no_hxc': False,         # Independent electron approximation
        'lrpa': False,           # Random Phase Approximation (no XC)
        'scissor': 0.0,          # Scissor shift in Rydberg units
        'charge_response': 0,    # Set to 1 to compute charge density response
        'pseudo_hermitian': True,  # Use pseudo-Hermitian Lanczos algorithm
        'd0psi_rs': False,       # Dipole computed in real space
        'lshift_d0psi': True     # Shift position operator for periodicity
    },
    'lr_post': {
        'omeg': 0.0,              # Transition energy in Rydberg units
        'epsil': 0.02,            # Broadening/damping term in Rydberg units
        'beta_gamma_z_prefix': 'pwscf',  # Prefix for beta gamma zeta coefficients
        'w_T_npol': 3,            # Number of polarization directions in previous calc
        'plot_type': 1            # Format for charge density response output
    }
}


mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run turbo_lanczos.x
qe_lanczosx(mof, lanczos_config)
