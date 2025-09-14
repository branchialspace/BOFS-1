# Run QuantumESPRESSO turbo_magnon.x

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
import numpy as np
from ase.io import read


def magnonx(
    structure_path,
    config
):
    """
    Run QE turbo_magnon.x calculation for magnon properties.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_magnonx_input(config, input_filename):
        """
        Write the QE turbo_magnon.x input file from config settings.
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
    calculation = "magnon"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['lr_input']['prefix'] = structure_name
    config['lr_input']['outdir'] = structure_name
    config['lr_control']['q1'] = config['q_vector'][0]
    config['lr_control']['q2'] = config['q_vector'][1]
    config['lr_control']['q3'] = config['q_vector'][2]
    # Write QE turbo_magnon.x input file
    write_magnonx_input(config, f"{run_name}.magi")
    # Subprocess run
    try:
        with open(f"{run_name}.mago", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.magi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("Magnon calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running magnon calculation: {cpe}")
        try:
            with open(f"{run_name}.mago", 'r') as f_out:
                print("\nMagnon Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


magnonx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_magnon.x'],
    'q_vector': [0.0, 0.0, 0.0],    # q-vector for magnon calculations in units of 2pi/a
    'lr_control': {
        'itermax': 500,             # Number of Lanczos iterations
        'pseudo_hermitian': True,   # Use pseudo-Hermitian Lanczos algorithm
        'approximation': 'TDDFT',   # Level of theory: 'TDDFT' or 'IPA'
        'ipol': 4,                  # Polarization: 1-3 for specific component, 4 for full tensor
    }
}


mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run turbo_magnon.x
qe_magnonx(mof, magnonx_config)
