# Ligand DFT and orbital analysis in ORCA
# Mayer, MBIS, NBO

import subprocess
import re
from ase import Atoms


def ligand_electron_analysis(ligand, charge=0, mult=1, method='B3LYP', basis_set='def2-SVP'):
    """

    Parameters:
    - ligand: ASE Atoms object representing the ligand molecule.
    - charge: Overall charge of the ligand.
    - mult: Spin multiplicity of the ligand.
    - method: Density functional method to use (e.g., 'B3LYP').
    - basis_set: Basis set to use (e.g., 'def2-SVP').

    Returns:
    - ligand_electron_analysis: Dictionary containing information about each atoms electron properties.
    """
    # Define the ORCA executable path
    orca_path = '/root/orca_6_0_1/orca'

    # Define filenames
    input_filename = 'orca_input.inp'
    output_filename = 'orca.out'

    # Write the ORCA input file
    write_orca_input(ligand, input_filename, charge, mult, method, basis_set)

    # Run ORCA
    try:
        with open(output_filename, 'w') as f_out:
            subprocess.run([orca_path, input_filename], stdout=f_out, stderr=subprocess.STDOUT, check=True)
        print("ORCA calculation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during ORCA calculation: {e}")
        return {}

    # Parse the energy
    energy = parse_energy(output_filename)
    print(f"Calculation completed. Energy: {energy}")

    # Parse the output file
    mayer_data = parse_mayer_data(output_filename)
    mbis_data = parse_mbis_data(output_filename)
    nbo_data = parse_nbo(output_filename)
    ligand_electron_analysis = process_orca_data(mayer_data, mbis_data, nbo_data)

    return ligand_electron_analysis

def write_orca_input(ligand, filename, charge, mult, method, basis_set):
    with open(filename, 'w') as f:
        # Write method and basis set
        f.write(f"! {method} {basis_set} MBIS NBO\n")
        f.write("%output\n")
        f.write("    Print[P_MBIS] 1\n")
        f.write("    Print[P_Mayer] 1\n")
        f.write("    Print[P_NBO] 1\n")
        f.write("    Print[P_NPA] 1\n")
        f.write("end\n\n")
        f.write("%method\n")
        f.write("    MAYER_BONDORDERTHRESH 0.05\n")
        f.write("end\n\n")
        # Write charge and multiplicity
        f.write(f"* xyz {charge} {mult}\n")
        # Write atomic coordinates
        for atom in ligand:
            symbol = atom.symbol
            x, y, z = atom.position
            f.write(f"  {symbol} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("*\n")

def parse_energy(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    energy = None
    for line in lines:
        if 'FINAL SINGLE POINT ENERGY' in line:
            parts = line.strip().split()
            if len(parts) >= 5:
                energy = float(parts[4])
                break
    return energy

def parse_mayer_data(filename):
    mayer_data = {'bond_orders': {}}
    with open(filename, 'r') as f:
        lines = f.readlines()

    bond_order_section = False
    bond_order_lines = []
    for i, line in enumerate(lines):
        if 'MAYER POPULATION ANALYSIS' in line:
            bond_order_section = True
            continue
        if bond_order_section and 'Mayer bond orders larger than' in line:
            # Start reading bond orders
            for j in range(i+1, len(lines)):
                line_content = lines[j].strip()
                if line_content == '':
                    break
                bond_order_lines.append(line_content)
            break

    # Combine bond order lines into a single string for regex parsing
    bond_order_text = ' '.join(bond_order_lines)

    # Now parse the bond order text
    bond_order_pattern = r'B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s*([\d\.]+)'
    matches = re.findall(bond_order_pattern, bond_order_text)
    for match in matches:
        atom1 = int(match[0])
        atom2 = int(match[1])
        order = float(match[2])
        mayer_data['bond_orders'][(atom1, atom2)] = order

    return mayer_data

def parse_mbis_data(filename):
    mbis_data = {}
    with open(filename, 'r') as f:
        lines = f.readlines()

    mbis_section = False
    valence_section = False
    for i, line in enumerate(lines):
        line = line.strip()
        if 'MBIS ANALYSIS' in line:
            mbis_section = True
            continue
        if mbis_section and 'ATOM     CHARGE    POPULATION     SPIN' in line:
            # Read total atomic charges and populations
            for j in range(i+1, len(lines)):
                line_content = lines[j].strip()
                if line_content == '':
                    break
                parts = line_content.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    atom_index = int(parts[0])
                    element = parts[1]
                    charge = float(parts[2])
                    population = float(parts[3])
                    if atom_index not in mbis_data:
                        mbis_data[atom_index] = {'element': element}
                    mbis_data[atom_index]['mbis_charge'] = charge
                    mbis_data[atom_index]['population'] = population
            continue

        if 'MBIS VALENCE-SHELL DATA:' in line:
            valence_section = True
            continue
        if valence_section and 'ATOM   POPULATION   WIDTH(A.U.)' in line:
            # Read valence shell data
            for j in range(i+1, len(lines)):
                line_content = lines[j].strip()
                if line_content == '':
                    valence_section = False
                    break
                parts = line_content.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    atom_index = int(parts[0])
                    element = parts[1]
                    valence_population = float(parts[2])
                    # width = float(parts[3])  # We can ignore width
                    if atom_index not in mbis_data:
                        mbis_data[atom_index] = {'element': element}
                    mbis_data[atom_index]['valence_population'] = valence_population
            continue
    return mbis_data

def parse_nbo(filename):
    nbo_data = {'npa_charges': {}, 'lone_pairs': {}}    
    with open(filename, 'r') as f:
        lines = f.readlines()

    npa_section = False
    lone_pair_section = False

    # Parse NPA charges
    for line in lines:
        if 'Summary of Natural Population Analysis:' in line:
            npa_section = True
            continue
        if npa_section:
            if '===' in line:
                npa_section = False
                continue
            parts = line.strip().split()
            if len(parts) >= 6 and parts[0].isalpha() and parts[1].isdigit():
                atom_index = int(parts[1]) - 1
                charge = float(parts[2])
                nbo_data['npa_charges'][atom_index] = charge

    # Parse lone pairs from the NBO summary
    lp_line_pattern = re.compile(r'^\s*(\d+\.)\s*LP\s*\(\s*(\d+)\s*\)\s+([A-Za-z]+)\s+(\d+)\s+([-\d\.]+)\s+([-\d\.]+)')
    for line in lines:
        if 'NATURAL BOND ORBITALS (Summary):' in line:
            lone_pair_section = True
            continue
        if lone_pair_section:
            if '------ non-Lewis' in line:
                break
            match = lp_line_pattern.match(line)
            if match:
                # Extract groups from regex, example:
                # group(1) = "9.", group(2) = LP number, group(3)=Element, group(4)=Atom index, group(5)=Occupancy, group(6)=Energy
                lp_number = int(match.group(2))
                element = match.group(3)
                atom_index = int(match.group(4)) - 1
                occupancy = float(match.group(5))
                energy = float(match.group(6))

                if atom_index not in nbo_data['lone_pairs']:
                    nbo_data['lone_pairs'][atom_index] = []
                
                nbo_data['lone_pairs'][atom_index].append({
                    'lp_number': lp_number,
                    'occupancy': occupancy,
                    'energy': energy
                })

    return nbo_data

def process_orca_data(mayer_data, mbis_data, nbo_data):
    """
    Return all data including NBO analysis.
    
    Parameters:
    - mayer_data: Dictionary containing Mayer bond orders
    - mbis_data: Dictionary containing MBIS charges and valence populations
    - nbo_data: Dictionary containing NPA charges and lone pair information
    
    Returns:
    - ligand_electron_analysis: Dictionary containing information about each atoms electron properties.
    """
    ligand_electron_analysis = {}

    for atom_index, mbis_info in mbis_data.items():
        element = mbis_info['element']
        mbis_charge = mbis_info.get('mbis_charge', None)
        valence_population = mbis_info.get('valence_population', None)

        # Get bonded electrons from Mayer bond orders
        bonded_electrons = sum(order for (a1, a2), order in mayer_data['bond_orders'].items() 
                             if atom_index in (a1, a2))

        # Get non-bonded electrons
        non_bonded_electrons = valence_population - bonded_electrons if valence_population is not None else None

        # Get NBO information
        npa_charge = nbo_data['npa_charges'].get(atom_index, None)
        lone_pairs = nbo_data['lone_pairs'].get(atom_index, [])
        
        # Calculate total lone pair electrons
        total_lp_electrons = sum(lp['occupancy'] for lp in lone_pairs)

        ligand_electron_analysis[atom_index] = {
            'element': element,
            'valence_population': valence_population,
            'mbis_charge': mbis_charge,
            'bonded_electrons': bonded_electrons,
            'non_bonded_electrons': non_bonded_electrons,
            'npa_charge': npa_charge,
            'nbo_lp': lone_pairs,
            'total_lp_electrons': total_lp_electrons
        }

    return ligand_electron_analysis
