# Run QuantumESPRESSO ph.x

import os
import subprocess
from subprocess import CalledProcessError
from pathlib import Path
from math import ceil, pi
import numpy as np
from ase.io import read


def phx(
    structure_path,
    config
):
    """
    Run QE ph.x calculation for phonon properties.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def qpoints(structure, q_spacing=0.25):
        """
        Given a desired q-point spacing q_spacing (in Å^-1),
        compute a suitable (nq1, nq2, nq3) Monkhorst–Pack grid for phonons.
        q_spacing : float
            Target spacing in reciprocal space, in Å^-1.
            For phonons, typically coarser than k-points.
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

    def write_phx_input(config, input_filename):
        """
        Write the QE ph.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # INPUTPH namelist
            f.write('&INPUTPH\n')
            for key, value in config['inputph'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')
            # Q-point specifications if not using ldisp or qplot
            if not config['inputph'].get('ldisp', False) and not config['inputph'].get('qplot', False):
                if 'xq' in config:
                    f.write(f"\n{config['xq'][0]:.8f} {config['xq'][1]:.8f} {config['xq'][2]:.8f}\n")
            # Q-points for qplot
            elif config['inputph'].get('qplot', False):
                if 'qpoints' in config:
                    f.write(f"\n{len(config['qpoints'])}\n")
                    for q_point in config['qpoints']:
                        f.write(f"{q_point[0]:.8f} {q_point[1]:.8f} {q_point[2]:.8f} {q_point[3]}\n")
            # Atom selection if nat_todo is specified
            if 'nat_todo' in config['inputph'] and config['inputph']['nat_todo'] > 0:
                if 'atoms' in config:
                    f.write('\n')
                    f.write(' '.join(str(atom) for atom in config['atoms']))
                    f.write('\n')
    
    # Args
    structure = read(structure_path)  # ASE Atoms object
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "ph"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['inputph']['prefix'] = structure_name
    config['inputph']['outdir'] = structure_name
    # Set q-points grid if ldisp is True
    if config['inputph'].get('ldisp', False):
        q_spacing = config.get('qpts_q_spacing', 0.25)
        nq1, nq2, nq3 = qpoints(structure, q_spacing)
        config['inputph']['nq1'] = nq1
        config['inputph']['nq2'] = nq2
        config['inputph']['nq3'] = nq3
    # Write QE ph.x input file
    write_phx_input(config, f"{run_name}.phi")
    # Subprocess run
    try:
        with open(f"{run_name}.pho", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.phi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("Phonon calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running phonon calculation: {cpe}")
        try:
            with open(f"{run_name}.pho", 'r') as f_out:
                print("\nPhonon Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


phx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/ph.x'],
    'xq': [0.0, 0.0, 0.0],         # q-point for non-ldisp calculations
    'qpts_q_spacing': 0.25,        # q-point spacing for automatic grid generation (coarser than k-points)
    'inputph': {
        'tr2_ph': 1.0e-14,         # Convergence threshold for phonons
        'ldisp': True,             # Run phonons on a grid of q-points
        'epsil': False,            # Calculate dielectric constant
        'trans': True,             # Calculate phonons
        'elop': True,              # Calculate electro-optic tensor
        'electron_phonon': '',     # electron-phonon coefficient method
        'fildyn': 'dynmat',        # Prefix for dynamical matrices
        'fildrho': 'drho',         # File for charge density response
        'fildvscf': 'dvscf',       # File for potential variation
        'max_seconds': 86400,      # Maximum allowed run time in seconds
        'asr': True,               # Apply Acoustic Sum Rule
        'search_sym': True,        # Enable mode symmetry analysis
    }
}

# ASE structure
mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run ph.x
qe_phx(mof, phx_config)
