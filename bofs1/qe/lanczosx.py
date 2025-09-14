# Run QuantumESPRESSO turbo_lanczos.x

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
from ase.io import read


def lanczosx(
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
    def write_lanczosx_input(config, input_filename):
        """
        Write the QE turbo_lanczos.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # input, control, post
            sections = ['lr_input', 'lr_control']
            if config['lr_control'].get('charge_response', 0) == 1:
                sections.append('lr_post')
            for section in sections:
                qe_section = section.upper()
                f.write(f'&{qe_section}\n')
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
    calculation = "turbo_lanczos"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['lr_input']['prefix'] = structure_name
    config['lr_input']['outdir'] = structure_name
    # Write QE turbo_lanczos.x input file
    write_lanczosx_input(config, f"{run_name}.tdfi")
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
