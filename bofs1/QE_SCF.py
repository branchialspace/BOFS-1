# QuantumESPRESSO Single-Point SCF with SOC

import os
import re
from pathlib import Path
from math import ceil, pi
import numpy as np
import subprocess
from subprocess import CalledProcessError
from ase import Atoms
from ase.io import read
from ase.data import atomic_masses, atomic_numbers
from ase.build import bulk
from mendeleev import element


def run_qe(
    structure,
    config
):
    """
    Run a QE calculation.
    
    structure : ase.Atoms
        An ASE Atoms object defining the system to calculate.
    config : dict
        Configuration dictionary containing required settings.
    """
    def pseudopotentials(structure, pseudo_directory):
        """
        Determine the appropriate pseudopotential file for each atomic species of the structure
        from a given directory. Tries to pick up 'fr' (fully relativistic) files
        if they are present; otherwise, picks the first found.
        """
        species = set(atom.symbol for atom in structure)
        pseudo_dict = {}
        pseudo_path = Path(pseudo_directory)
        for symbol in species:
            candidates = []
            for pp_file in pseudo_path.glob('*.[uU][pP][fF]'):
                pp_stem = pp_file.stem.split('_')[0].split('.')[0]
                if pp_stem.lower() == symbol.lower():
                    candidates.append(pp_file)
            if not candidates:
                raise FileNotFoundError(f"No pseudopotential found for {symbol}")
            # Prefer 'fr' if available
            fr_files = [c for c in candidates if 'fr' in c.stem.lower()]
            selected = max(fr_files, key=lambda x: len(x.stem)) if fr_files else candidates[0]
            pseudo_dict[symbol] = selected.name

        return pseudo_dict

    def pseudo_cutoffs(pseudo_dict, pseudo_directory, wfn_scalar=1.15, rho_scalar=1.15):
        """
        Parse pseudopotential files for wavefunction and density cutoffs, 
        select the largest and scale.

        wfn_scalar : float
            Scaling factor for wavefunction cutoff.
        rho_scalar : float
            Scaling factor for density cutoff.

        Returns 
        ecutwfc, ecutrho : tuple
            containing the largest scaled and
            rounded cutoff values for wavefunctions and density.
        """
        all_wfc = []
        all_rho = []
        for _, pp_filename in pseudo_dict.items():
            pp_path = Path(pseudo_directory) / pp_filename
            with open(pp_path, 'r') as f:
                content = f.read()
            # Regex to grab wfc_cutoff, rho_cutoff
            wfc_match = re.search(r'wfc_cutoff\s*=\s*"?\s*(?P<val>[\d.+Ee-]+)\s*"?', content)
            rho_match = re.search(r'rho_cutoff\s*=\s*"?\s*(?P<val>[\d.+Ee-]+)\s*"?', content)
            # Parse wfc_cutoff if found
            if wfc_match:
                raw_wfc = wfc_match.group('val').strip()
                # Normalize exponentials
                raw_wfc = raw_wfc.replace('E+', 'E').replace('e+', 'e')
                wfc_val = float(raw_wfc)
                all_wfc.append(wfc_val)
            # Parse rho_cutoff if found
            if rho_match:
                raw_rho = rho_match.group('val').strip()
                # Normalize exponentials
                raw_rho = raw_rho.replace('E+', 'E').replace('e+', 'e')
                rho_val = float(raw_rho)
                all_rho.append(rho_val)
        # Decide final cutoffs based on presence of wfc_cutoff
        if len(all_wfc) == 0:
            # That means each file must have had only 'rho_cutoff', which is wfc.
            # We'll use the maximum of those as ecutwfc, and ecutrho = 4 * ecutwfc.
            if not all_rho:
                # No cutoffs found at all
                raise ValueError("No wfc_cutoff or rho_cutoff found in any pseudopotential file.")
            max_wfc = max(all_rho)
            max_rho = 4.0 * max_wfc
        else:
            # We have proper wfc_cutoff in each file (and presumably also rho_cutoff).
            max_wfc = max(all_wfc)
            max_rho = max(all_rho)
        # Multipliers to be safe, then round up
        ecutwfc = ceil(max_wfc * wfn_scalar)
        ecutrho = ceil(max_rho * rho_scalar)

        return ecutwfc, ecutrho

    def kpoints(structure, k_spacing=0.13, shift=(1,1,1)):
        """
        Given a desired k-point spacing k_spacing (in Å^-1),
        compute a suitable (n1, n2, n3) Monkhorst–Pack grid for the structure.

        k_spacing : float
            Target spacing in reciprocal space, in Å^-1.
        shift : tuple of int
            The Monkhorst-Pack shift for each direction, either (0,0,0) or (1,1,1).

        Returns
        (n1, n2, n3, s1, s2, s3) : tuple of ints
            The grid subdivisions (n1, n2, n3) and the shift (s1, s2, s3).
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
        n1 = max(1, ceil(b1_len / k_spacing))
        n2 = max(1, ceil(b2_len / k_spacing))
        n3 = max(1, ceil(b3_len / k_spacing))

        # Unpack the shift
        s1, s2, s3 = shift

        return (n1, n2, n3, s1, s2, s3)

    def nbnd(structure, nbnd_scalar = 2):
        """
        Number of electronic states (bands) to be calculated,
        scaled total valence count from Mendeleev.

        nbnd_scalar: float
            Scaling factor for bands.
        
        Returns
        nbnd : int
            Number of electronic states (bands).
        """
        # Sum valence electrons for each atom
        total_valence = 0
        for atom in structure:
            symbol = atom.symbol
            elem = element(symbol)
            valence = elem.nvalence()
            total_valence += valence
        # Calculate recommended nbnd with safety factor
        # Typically we use max(int(total_valence * 0.7), int(total_valence * 0.5) + 4) for a metal
        nbnd = (int(total_valence * nbnd_scalar))

        return nbnd

    def write_espresso_input(
        structure,
        config,
        pseudopotentials,
        kpoints,
        input_filename
    ):
        """
        Write the QE input file from an ASE structure, config,
        and selected pseudopotentials.
        """
        with open(input_filename, 'w') as f:
            # Control section
            f.write('&control\n')
            for key, value in config['control'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')

            # System section
            f.write('&system\n')
            for key, value in config['system'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')

            # Electrons section
            f.write('&electrons\n')
            for key, value in config['electrons'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')

            # Atomic species
            f.write('\nATOMIC_SPECIES\n')
            unique_symbols = set(structure.get_chemical_symbols())
            for symbol in unique_symbols:
                mass = atomic_masses[atomic_numbers[symbol]]
                f.write(f"  {symbol.title()} {mass:.4f} {pseudopotentials[symbol]}\n")

            # Atomic positions
            f.write('\nATOMIC_POSITIONS angstrom\n')
            abs_pos = structure.get_positions()
            for atom, pos in zip(structure, abs_pos):
                f.write(f"  {atom.symbol.title()} "
                        f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")

            # K-points grid
            f.write('\nK_POINTS automatic\n')
            f.write(f"  {kpoints[0]} {kpoints[1]} {kpoints[2]} {kpoints[3]} {kpoints[4]} {kpoints[5]}\n")

            # Cell parameters in angstrom
            f.write('\nCELL_PARAMETERS angstrom\n')
            for vec in structure.cell:
                f.write(f"  {vec[0]:.10f} {vec[1]:.10f} {vec[2]:.10f}\n")

    # Unpack top-level config items
    command = config['command']
    pseudo_dir = config['pseudo_dir']
    input_filename = config['input_filename']
    output_filename = config['output_filename']
    outdir = config['outdir']
    os.makedirs(outdir, exist_ok=True)
    # Set nat and ntyp based on the current structure
    config['system']['nat'] = len(structure)
    config['system']['ntyp'] = len(set(structure.get_chemical_symbols()))
    # Set Pseudos
    pseudopotentials = pseudopotentials(structure, pseudo_dir)
    # Set Cutoffs
    wfn_scalar = config['wfn_scalar']
    rho_scalar = config['rho_scalar']
    ecutwfc, ecutrho = pseudo_cutoffs(pseudopotentials, pseudo_dir)
    config['system']['ecutwfc'] = ecutwfc
    config['system']['ecutrho'] = ecutrho
    # Set k-points
    k_spacing = config['kpts_k_spacing']
    shift = config['kpts_shift']
    kpoints = kpoints(structure, k_spacing, shift)
    # Set nbnd
    nbnd_scalar = config['nbnd_scalar']
    config['system']['nbnd'] = nbnd(structure, nbnd_scalar)
    # write QE input file
    write_espresso_input(structure, config, pseudopotentials, kpoints, input_filename)
    # Subprocess run
    try:
        with open(output_filename, 'w') as f_out:
            command_list = config['command'] + ['-in', input_filename]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True
            )
        print("QE calculation completed successfully.")

    except CalledProcessError as cpe:
        print(f"Error running QE: {cpe}")
        try:
            with open(output_filename, 'r') as f_out:
                print("\nQE Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")

    except Exception as e:
        print(f"Unexpected error: {e}")


config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/pw.x'],
    'pseudo_dir': '/content/ONCVPSP/abinit/',
    'input_filename': 'espresso.pwi',
    'output_filename': 'espresso.pwo',
    'outdir': '/content/tmp',
    'wfn_scalar': 1.15,
    'rho_scalar': 1.15,
    'kpts_k_spacing': 0.13,
    'kpts_shift': (1,1,1),
    'nbnd_scalar': 2,
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': '/content/ONCVPSP/abinit/',
        'outdir': '/content/tmp',
        'prefix': 'TIRDOO_full',
        'disk_io': 'medium',
        'wf_collect': True,
        'tprnfor': True,
        'tstress': True,
    },
    'system': {
        'ibrav': 0,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,
        'noncolin': True,
        'lspinorb': True
    },
    'electrons': {
        'conv_thr': 1.0e-6,
        'mixing_beta': 0.3,
        'electron_maxstep': 300,
    }
}

# ASE structure
# bi = bulk('Bi', 'rhombohedral', a=4.75, c=12.36, orthorhombic=False)
mof = read("/content/mofs/TIRDOO_full.cif")
# Run QE
run_qe(mof, config)
