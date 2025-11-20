# Firtst version of the Hubbard card function for pwx with blockwise orbital selection
def hubbard_atoms(structure, pseudo_dict, pseudo_directory, initial_u_value=0.1, n_manifolds=1):
    """
    Identify atoms needing Hubbard U + V corrections, extract valence orbitals
    from pseudopotential files, and prioritize manifolds based on orbital blocks.
    For initial run only, U + V pairs and values will be inferred by hp.x for subsequent runs.
    structure : ASE Atoms object
        The atomic structure
    pseudo_dict : dict
        Dictionary mapping atom symbols to pseudopotential filenames
    pseudo_directory : str
        Path to the directory containing pseudopotential files
    initial_u_value : float
        Initial U value to assign (will be refined by hp.x)
    n_manifolds : int
        Number of orbitals/ manifolds per species
    hubbard_card : list
        List with manifold information formatted for the hubbard card
    """
    # Species known to never require Hubbard corrections
    non_correlated_species = {'H', 'He', 'Li', 'Be', 'B', 'F', 'Ne',
                              'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar',
                              'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                              'Rb', 'Sr', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
                              'Cs', 'Ba', 'Hg', 'Tl', 'Po', 'At', 'Rn'}
    # Get unique atom types in structure
    atom_types = sorted(set(structure.get_chemical_symbols()))
    # Identify Hubbard candidates
    hubbard_candidates = [sym for sym in atom_types if sym not in non_correlated_species]
    # Write Hubbard card
    hubbard_card = []
    # Extract orbital labels from pseudopotential file
    for symbol in hubbard_candidates:
        pp_filename = pseudo_dict[symbol]
        pp_path = Path(pseudo_directory) / pp_filename
        with open(pp_path, 'r') as f:
            content = f.read()
        # Find all PP_CHI entries with label
        pattern = r'<PP_CHI\.\d+.*?label\s*=\s*"([^"]+)"'
        matches = re.findall(pattern, content, re.DOTALL)
        orbital_info = []
        # Parse orbital labels to extract (n) principal quantum number and (l) angular momentum 
        for label in matches:
            n = int(''.join(filter(str.isdigit, label))) if any(c.isdigit() for c in label) else None
            l_type = next((char for char in label if char in 'spdf'), None)
            if n is not None and l_type is not None:
                orbital_info.append({'label': label, 'n': n, 'l_type': l_type})
        orbital_info = list({orb['label']: orb for orb in orbital_info}.values())
        # Use Mendeleev to determine element type and electron orbital block
        elem = element(symbol)
        block = elem.block
        # Sort orbitals by type
        s_orbitals = sorted([o for o in orbital_info if o['l_type'] == 's'], key=lambda o: -o['n'])
        p_orbitals = sorted([o for o in orbital_info if o['l_type'] == 'p'], key=lambda o: -o['n'])
        d_orbitals = sorted([o for o in orbital_info if o['l_type'] == 'd'], key=lambda o: -o['n'])
        f_orbitals = sorted([o for o in orbital_info if o['l_type'] == 'f'], key=lambda o: -o['n'])
        # Prioritize orbitals by block of species, special case for oxygen
        if symbol == 'O':
            o_2p = [o for o in p_orbitals if o['n'] == 2]
            other_p = [o for o in p_orbitals if o['n'] != 2]
            prioritized_orbitals = o_2p + other_p + s_orbitals + d_orbitals + f_orbitals
        elif block == 'd':
            prioritized_orbitals = d_orbitals + p_orbitals + s_orbitals + f_orbitals
        elif block == 'p':
            prioritized_orbitals = p_orbitals + s_orbitals + d_orbitals + f_orbitals
        elif block == 'f':
            prioritized_orbitals = f_orbitals + d_orbitals + p_orbitals + s_orbitals
        elif block == 's':
            prioritized_orbitals = s_orbitals + p_orbitals + d_orbitals + f_orbitals
        top_manifolds = prioritized_orbitals[:min(n_manifolds, len(prioritized_orbitals))]
        # Format Hubbard card
        if top_manifolds:
            # First manifold
            hubbard_card.append(f"U {symbol}-{top_manifolds[0]['label']} {initial_u_value:.1f}")
            # Combine second and third manifolds if they exist
            if len(top_manifolds) > 1:
                combined_labels = '-'.join(orb['label'] for orb in top_manifolds[1:])
                hubbard_card.append(f"U {symbol}-{combined_labels} {initial_u_value:.1f}")

    return hubbard_card
