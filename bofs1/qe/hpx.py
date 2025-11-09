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
    def qpoints(structure, q_spacing=0.065):
        """
        Given a desired q-point spacing q_spacing (in Å^-1),
        compute a suitable (nq1, nq2, nq3) Monkhorst–Pack grid for Hubbard parameters.
        q_spacing : float
            Target spacing in reciprocal space, in Å^-1.
            For Hubbard parameters, typically denser than k-points.
        Returns
        (nq1, nq2, nq3) : tuple of ints
            The grid subdivisions.
        """
        # Extract real-space lattice vectors
        cell = structure.get_cell()  # 3x3 array
        a1, a2, a3 = [np.array(vec) for vec in cell]
        # Compute real-space volume
        volume = np.dot(a1, np.cross(a2, a3))
        # Compute reciprocal lattice vectors b1, b2, b3
        # b1 = 2π * (a2 × a3) / (a1 · (a2 × a3)), etc.
        b1 = 2 * pi * np.cross(a2, a3) / volume
        b2 = 2 * pi * np.cross(a3, a1) / volume
        b3 = 2 * pi * np.cross(a1, a2) / volume
        # Compute magnitudes of reciprocal vectors
        b1_len = np.linalg.norm(b1)
        b2_len = np.linalg.norm(b2)
        b3_len = np.linalg.norm(b3)
        # Determine the number of divisions along each direction
        # Small reciprocal lattice vectors (in Å⁻¹) indicate large unit cell dims
        dim_threshold = 0.05  # threshold in Å⁻¹, corresponds to ~125Å real-space dimension
        n1 = max(1, ceil(b1_len / q_spacing)) if b1_len > dim_threshold else 1
        n2 = max(1, ceil(b2_len / q_spacing)) if b2_len > dim_threshold else 1
        n3 = max(1, ceil(b3_len / q_spacing)) if b3_len > dim_threshold else 1

        return (n1, n2, n3)

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
            f.write('/')
    # Args
    structure = read(structure_path)  # ASE Atoms object
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "hp"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['inputhp']['prefix'] = structure_name
    config['inputhp']['outdir'] = structure_name
    # Set q-points
    q_spacing = config['qpts_q_spacing']
    nq1, nq2, nq3 = qpoints(structure, q_spacing)
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
