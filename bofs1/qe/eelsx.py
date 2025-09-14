# Run QuantumESPRESSO turbo_eels.x

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
import numpy as np
from ase.io import read


def eelsx(
    structure_path,
    config
):
    """
    Run QE turbo_eels.x calculation for electron energy loss spectra.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_eelsx_input(config, input_filename):
        """
        Write the QE turbo_eels.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # input, control
            for section in ['lr_input', 'lr_control']:
                section_name = section.upper()
                f.write(f'&{section_name}\n')
                for key, value in config[section].items():
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
    calculation = "eels"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['lr_input']['prefix'] = structure_name
    config['lr_input']['outdir'] = structure_name
    config['lr_control']['q1'] = config['q_vector'][0]
    config['lr_control']['q2'] = config['q_vector'][1]
    config['lr_control']['q3'] = config['q_vector'][2]
    # Write QE turbo_eels.x input file
    write_eelsx_input(config, f"{run_name}.eelsi")
    # Subprocess run
    try:
        with open(f"{run_name}.eelso", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.eelsi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("EELS calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running EELS calculation: {cpe}")
        try:
            with open(f"{run_name}.eelso", 'r') as f_out:
                print("\nEELS Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


eelsx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_eels.x'],
    'q_vector': [0.0, 0.0, 0.0],    # q-vector for EELS calculations in units of 2pi/a
    'lr_control': {
        'approximation': 'TDDFT',   # Level of theory: 'TDDFT', 'IPA', or 'RPA_with_CLFE'
        'calculator': 'lanczos',    # Algorithm: 'lanczos' or 'sternheimer'
        'itermax': 500,             # Number of iterations (Lanczos or Sternheimer)
        'pseudo_hermitian': True,   # Use pseudo-Hermitian Lanczos algorithm
        'ethr_nscf': 1.0e-11        # Threshold for convergence of eigenvalues
        # The following parameters are used when calculator = 'sternheimer'
        # 'alpha_mix': 0.7,         # Mixing parameter for SCF potential
        # 'epsil': 0.02,            # Broadening/damping term in Rydberg units
        # 'units': 1,               # Units: 0=Rydbergs, 1=Electron volts
        # 'start': 0.0,             # Starting frequency
        # 'end': 2.5,               # Ending frequency
        # 'increment': 0.001        # Frequency step
    }
}


mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run turbo_eels.x
qe_eelsx(mof, eelsx_config)
