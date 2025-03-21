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


def qe_hpx(
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
    def qpoints(structure, q_spacing=0.25):
        """
        Given a desired q-point spacing q_spacing (in Å^-1),
        compute a suitable (nq1, nq2, nq3) Monkhorst–Pack grid for Hubbard parameters.
        q_spacing : float
            Target spacing in reciprocal space, in Å^-1.
            For Hubbard parameters, typically coarser than k-points (0.2-0.3 is common).
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
    
    def hubbard_atoms(structure):
        """
        Identify atoms needing Hubbard U/V corrections, excluding species
        that are definitively non-correlated.
        Returns
        hubbard_atoms : dict
            skip_type : List of booleans for each atom type
            equiv_type : List of equivalence indices for each atom type
            perturb_only_atom : List of booleans for each atom type
            hubbard_candidates : List of (index, symbol) tuples for atoms needing correction
        """
        # Species known to never require Hubbard corrections
        non_correlated_species = {'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne',
                                  'Si', 'P', 'S', 'Cl', 'Ar', 'Ge', 'As', 'Se', 'Br', 'Kr',
                                  'I', 'Xe', 'Rn'}
        # Get unique atom types in structure
        atom_types = sorted(set(structure.get_chemical_symbols()))
        n_types = len(atom_types)
        # Initialize parameters
        skip_type = [symbol in non_correlated_species for symbol in atom_types]
        equiv_type = [0] * n_types
        perturb_only_atom = [False] * n_types
        # Only atoms not skipped are considered Hubbard candidates
        hubbard_candidates = [(i, symbol) for i, symbol in enumerate(atom_types) if not skip_type[i]]
        hubbard_atoms = {
            'skip_type': skip_type,
            'equiv_type': equiv_type,
            'perturb_only_atom': perturb_only_atom,
            'hubbard_candidates': hubbard_candidates
        }

        return hubbard_atoms

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
            # Hubbard atoms
            for i, (skip, equiv, perturb) in enumerate(zip(config['skip_type'], config['equiv_type'], config['perturb_only_atom']), 1):
                f.write(f"  skip_type({i}) = {'.true.' if skip else '.false.'}\n")
                f.write(f"  equiv_type({i}) = {equiv}\n")
                f.write(f"  perturb_only_atom({i}) = {'.true.' if perturb else '.false.'}\n")

            f.write('/\n')
        
    # Args
    structure = read(structure_path)  # ASE Atoms object
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "hp"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['inputhp']['prefix'] = structure_name
    config['inputhp']['outdir'] = structure_name
    os.makedirs(structure_name, exist_ok=True)
    # Set q-points
    q_spacing = config.get('qpts_q_spacing', 0.25)
    nq1, nq2, nq3 = qpoints(structure, q_spacing)
    config['inputhp']['nq1'] = nq1
    config['inputhp']['nq2'] = nq2
    config['inputhp']['nq3'] = nq3
    # Set Hubbard atoms
    hubbard_params = hubbard_atoms(structure)
    config['skip_type'] = hubbard_params['skip_type']
    config['equiv_type'] = hubbard_params['equiv_type']
    config['perturb_only_atom'] = hubbard_params['perturb_only_atom']
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


hpx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/hp.x'],
    'qpts_q_spacing': 0.25,        # q-point spacing
    'inputhp': {
        'iverbosity': 2,           # Verbosity level
        'find_atpert': 4,          # Method for determining which atoms to perturb (4=Perturb all Hubbard atoms)
        'docc_thr': 5.0e-5,        # Threshold for comparing unperturbed occupations
        'skip_equivalence_q': True, # If True, use full q-point grid; if False, use symmetry to reduce q-points
        'conv_thr_chi': 1.0e-5,    # Convergence threshold for response function chi
        'thresh_init': 1.0e-14,    # Initial threshold for linear system solution
        'ethr_nscf': 1.0e-11,      # Threshold for NSCF eigenvalue convergence
        'niter_max': 100,          # Maximum number of iterations
        'alpha_mix': 0.3,          # Mixing parameter for updating SCF potential
        'nmix': 4,                 # Number of iterations used in Broyden potential mixing
        'compute_hp': False,       # If True, post-process previously calculated results
        'determine_num_pert_only': False, # If True, determine only the number of perturbations and exit
        'determine_q_mesh_only': False,   # If True, determine only the q mesh and exit
        'sum_pertq': False,        # If True, collect response matrices for all q points
        'max_seconds': 86400,      # Maximum allowed run time in seconds
        'num_neigh': 6,            # Number of nearest neighbors for Hubbard V parameters
        'lmin': 2,                 # Minimum orbital quantum number for Hubbard V parameters
        'rmax': 100.0,             # Maximum distance to search for neighbors
        'dist_thr': 6.0e-4,        # Threshold for comparing inter-atomic distances
    }
}

mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run hp.x
qe_hpx(mof, hpx_config)
