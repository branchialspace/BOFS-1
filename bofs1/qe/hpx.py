# Run QuantumESPRESSO hp.x

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
from math import ceil, pi
import numpy as np
from ase.io import read
from ase.data import chemical_symbols
from mendeleev import element


def hpx(
    structure_path,
    config
):
    """
    Run QE hp.x calculation for Hubbard parameters.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def qpoints(prefix):
        """
        Read the SCF k-point mesh from PREFIX.pwi.
        Expects block:
        K_POINTS
          nk1 nk2 nk3 1 1 1
        """
        pwi_file = f"{prefix}.pwi"
        n1 = n2 = n3 = None
        with open(pwi_file, "r") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().lower().startswith("k_points"):
                parts = lines[i+1].split()
                n1, n2, n3 = map(int, parts[:3])
                break
    
        return n1, n2, n3

    def write_hpx_input(config, input_filename):
        """
        Write the QE hp.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # INPUTHP namelist
            f.write('&INPUTHP\n')
            for key, value in config['inputhp'].items():
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
    calculation = "hp"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['inputhp']['prefix'] = structure_name
    config['inputhp']['outdir'] = structure_name
    # Set q-points
    prefix = config['inputhp']['prefix']
    nq1, nq2, nq3 = qpoints(prefix)
    config['inputhp']['nq1'] = nq1
    config['inputhp']['nq2'] = nq2
    config['inputhp']['nq3'] = nq3
    # Write QE hp.x input file
    write_hpx_input(config, f"{run_name}.hpi")
    # Subprocess run
    try:
        with open(f"{run_name}.hpo", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.hpi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("Hubbard parameters calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running Hubbard parameters calculation: {cpe}")
        try:
            with open(f"{run_name}.hpo", 'r') as f_out:
                print("\nHubbard Parameters Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
