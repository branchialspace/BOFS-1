# Run QuantumESPRESSO Plane-Wave Self-Consistent Field

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


def qe_pwx(
    structure_path,
    config
):
    """
    Run QE pw.x PWscf calculation.
    structure_path : string
        Path of a cif unit cell defining the system to calculate.
    config : dict
        Configuration dictionary containing required settings.
    """
    def pseudopotentials(structure, pseudo_directory):
        """
        Determine the appropriate pseudopotential file for each atomic species of the structure
        from a given directory. Tries to pick up 'fr' (fully relativistic) files
        if they are present; otherwise, picks the first found.
        Returns
        pseudo_dict : dict
            Dictionary of atomic species and associated pseudopotential files.
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

    def charge(structure_path, structure_name):
        """
        Determine the total charge of the system from the structure name or (if '_charged'
        is present) by parsing the CIF file's `_chemical_formula_moiety` key.
        Returns : float
            The system's total charge.
        """
        charge = 0.0
        # If "_charged" is in structure_name
        if "_charged" in structure_name.lower():
            with open(structure_path, 'r') as f:
                formula_line = None
                for line in f:
                    if '_chemical_formula_moiety' in line:
                        match = re.search(r"_chemical_formula_moiety\s+'([^']+)'", line)
                        if match:
                            formula_line = match.group(1)
                            break
            if formula_line:
                match_charge = re.search(r'(\d+)(\+|\-)\)n', formula_line)
                if match_charge:
                    digits = match_charge.group(1)
                    sign_symbol = match_charge.group(2)
                    if sign_symbol == '+':
                        charge = float(digits)
                    else:
                        charge = -float(digits)
        else:
            # If structure_name ends with "_px" or "_nx" (x is integer)
            match_sign_digits = re.search(r'_([pn])(\d+)$', structure_name)
            if match_sign_digits:
                sign = match_sign_digits.group(1)
                digits = match_sign_digits.group(2)
                if sign == 'p':
                    charge = float(digits)
                else:
                    charge = -float(digits)

        return charge

    def hubbard_atoms(structure, pseudo_dict, pseudo_directory, initial_u_value=0.1):
        """
        Identify atoms needing Hubbard U+V corrections, extract valence orbitals
        from pseudopotential files, and prioritize manifolds based on orbital energy proximity to fermi level.
        For initial run only, U + V pairs and values will be inferred by hp.x for subsequent runs.
        structure : ASE Atoms object
            The atomic structure
        pseudo_dict : dict
            Dictionary mapping atom symbols to pseudopotential filenames
        pseudo_directory : str
            Path to the directory containing pseudopotential files
        initial_u_value : float
            Initial U value to assign (will be refined by hp.x)
        Returns
        hubbard_data : dict
            Dictionary with manifold information and formatted hubbard card
        """
        # Species known to never require Hubbard corrections
        non_correlated_species = {'H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne',
                                'Si', 'P', 'S', 'Cl', 'Ar', 'Ge', 'As', 'Se', 'Br', 'Kr',
                                'I', 'Xe', 'Rn'}
        # Get unique atom types in structure
        atom_types = sorted(set(structure.get_chemical_symbols()))
        # Identify Hubbard candidates
        hubbard_candidates = [symbol for symbol in atom_types if symbol not in non_correlated_species]
        # Store manifold information and formatted hubbard card
        hubbard_manifolds = {}
        hubbard_values = {}
        hubbard_card = []
        for symbol in hubbard_candidates:
            pp_filename = pseudo_dict[symbol]
            pp_path = Path(pseudo_directory) / pp_filename
            # Extract orbital information from pseudopotential file
            orbital_info = []
            with open(pp_path, 'r') as f:
                content = f.read()
            # Find all PP_CHI entries with label, l value, occupation, pseudo_energy
            pattern = r'<PP_CHI\.\d+.*?label\s*=\s*"([^"]+)".*?l\s*=\s*"(\d+)".*?occupation\s*=\s*"([^"]+)"(?:.*?pseudo_energy\s*=\s*"([^"]+)")?'
            matches = re.findall(pattern, content, re.DOTALL)
            for match in matches:
                label, l_value, occupation, energy = match  # energy is None in dalcorso paw pseudos
                # If energy is not available, get it from Mendeleev
                if not (energy and energy.strip()):
                    # Use l_value which is already in spectroscopic notation (0=s, 1=p, 2=d, 3=f)
                    l = int(l_value)
                    # Extract principal quantum number from label
                    n = None
                    if label and any(c.isdigit() for c in label):
                        # Extract all digits from the beginning of the label
                        digits = ""
                        for c in label:
                            if c.isdigit():
                                digits += c
                            else:
                                break
                        if digits:
                            n = int(digits)
                    # Get energy from Mendeleev
                    elem = element(symbol)
                    # Find the appropriate orbital energy
                    for orbital in elem.orbitals:
                        if orbital.l == l and (n is None or orbital.n == n):
                            energy = orbital.energy
                            break
                    else:
                        energy = 0.0  # Fallback if no matching orbital found
                else:
                    energy = float(energy.strip())
                orbital_info.append({
                    'label': label.lower(),
                    'l_value': int(l_value),
                    'occupation': float(occupation.strip()),
                    'energy': energy})
            # For duplicate orbitals, select instance with energy closest to Fermi level (assumed to be at zero)
            orbital_info = [sorted([o for o in orbital_info if o['label'] == label], key=lambda x: abs(x['energy']))[0] for label in {o['label'] for o in orbital_info}]
            # Sort orbitals by priority: proximity to Fermi level (assumed to be at zero energy)
            sorted_orbitals = sorted(
                orbital_info,
                key=lambda x: abs(x['energy']),  # Proximity to Fermi level (assumed to be at zero)
                reverse=False)  # Smaller absolute difference first
            # Get top manifold
            top_manifolds = sorted_orbitals[:min(1, len(sorted_orbitals))]
            hubbard_manifolds[symbol] = [orbital['label'] for orbital in top_manifolds]
            for orbital in top_manifolds:
                hubbard_values[(symbol, orbital['label'])] = initial_u_value
            # Format the Hubbard card entries
            if top_manifolds:
                # First manifold
                first_orbital = top_manifolds[0]
                hubbard_card.append(f"U {symbol}-{first_orbital['label']} {initial_u_value:.1f}")
                # Combine second and third manifolds if they exist
                if len(top_manifolds) > 1:
                    combined_labels = '-'.join([orbital['label'] for orbital in top_manifolds[1:]])
                    hubbard_card.append(f"U {symbol}-{combined_labels} {initial_u_value:.1f}")
        hubbard_data = {
            'hubbard_manifolds': hubbard_manifolds,
            'hubbard_values': hubbard_values,
            'hubbard_card': hubbard_card}
    
        return hubbard_data

    def write_pwx_input(
        structure,
        config,
        pseudopotentials,
        kpoints,
        input_filename
    ):
        """
        Write the QE PWscf input file from an ASE structure, config,
        and selected pseudopotentials.
        """
        with open(input_filename, 'w') as f:
            # Control, System, Electrons
            for section in ['control', 'system', 'electrons']:
                f.write(f'&{section}\n')
                for key, value in config[section].items():
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
            # Hubbard U+V corrections
            if hubbard_data and hubbard_data.get('hubbard_card'):
                f.write('\nHUBBARD ortho-atomic\n')
                for line in hubbard_data['hubbard_card']:
                    f.write(f"{line}\n")

    # Args
    structure = read(structure_path) # ASE Atoms object
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = config['control']['calculation']
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    pseudo_dir = config['control']['pseudo_dir']
    config['control']['prefix'] = structure_name
    config['control']['outdir'] = structure_name
    os.makedirs(structure_name, exist_ok=True)
    # Set nat and ntyp
    config['system']['nat'] = len(structure)
    config['system']['ntyp'] = len(set(structure.get_chemical_symbols()))
    # Set Pseudos
    pseudopotentials = pseudopotentials(structure, pseudo_dir)
    # Set Cutoffs
    wfn_scalar = config['wfn_scalar']
    rho_scalar = config['rho_scalar']
    ecutwfc, ecutrho = pseudo_cutoffs(pseudopotentials, pseudo_dir, wfn_scalar, rho_scalar)
    config['system']['ecutwfc'] = ecutwfc
    config['system']['ecutrho'] = ecutrho
    # Set k-points
    k_spacing = config['kpts_k_spacing']
    shift = config['kpts_shift']
    kpoints = kpoints(structure, k_spacing, shift)
    # Set nbnd
    nbnd_scalar = config['nbnd_scalar']
    config['system']['nbnd'] = nbnd(structure, nbnd_scalar)
    # Set total charge
    charge = charge(structure_path, structure_name)
    config['system']['tot_charge'] = charge
    # Set Hubbard corrections
    hubbard_data = hubbard_atoms(structure, pseudopotentials, pseudo_dir)
    # Write QE PWscf input file
    write_pwx_input(structure, config, pseudopotentials, kpoints, f"{run_name}.pwi")
    # Subprocess run
    try:
        with open(f"{run_name}.pwo", 'w') as f_out:
            command_list = config['command'] + ['-in', f"{run_name}.pwi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("QE calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running QE: {cpe}")
        try:
            with open(f"{run_name}.pwo", 'r') as f_out:
                print("\nQE Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/pw.x'],
    'wfn_scalar': 1.15,
    'rho_scalar': 1.15,
    'kpts_k_spacing': 0.13, # scf: 0.13    nscf: 0.09
    'kpts_shift': (1,1,1),
    'nbnd_scalar': 2.5,
    'control': {
        'calculation': 'scf', # scf     nscf     bands
        'restart_mode': 'from_scratch',
        'pseudo_dir': '/content/pslibrary/rel-pbe/PSEUDOPOTENTIALS', # /content/ONCVPseudoPack/Abinit_v0.4/UPF/PBEsol   /content/pslibrary/rel-pbe/PSEUDOPOTENTIALS
        'disk_io': 'medium',
        'verbosity': 'high',
        'wf_collect': True,
        'tprnfor': True,
        'tstress': True
    },
    'system': {
        'ibrav': 0,
        'occupations': 'smearing',
        'smearing': 'fermi-dirac', # gaussian     marzari-vanderbilt     fermi-dirac
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
mof = "/content/mofs/SIWZOO_full_n2.cif" # /content/mofs/SIWZOO_full_n2.cif      /content/mofs/TIRDOO_full.cif
# Run PWscf
qe_pwx(mof, config)
