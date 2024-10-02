# Mayer and Minimal Basis Iterative Stockholder (MBIS) electron-donor analysis in ORCA

import os
import subprocess
from ase import Atoms
from ase.data import chemical_symbols

def electron_donors(ligand, charge=0, mult=1, method='B3LYP', basis_set='def2-SVP'):
    # Define the ORCA executable path
    orca_path = '/root/orca_6_0_0/orca'

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
    donor_atoms = find_electron_donors(mayer_data, mbis_data)

    return donor_atoms

def write_orca_input(ligand, filename, charge, mult, method, basis_set):
    with open(filename, 'w') as f:
        # Write method and basis set
        f.write(f"! {method} {basis_set} MBIS\n")
        f.write("%output\n")
        f.write("    Print[P_MBIS] 1\n")
        f.write("    Print[P_Mayer] 1\n")
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
                    continue
                if '****' in line_content or 'Population analysis' in line_content:
                    # Reached the end of bond orders
                    break
                bond_order_lines.append(line_content)
            break

    # Combine bond order lines into a single string for regex parsing
    bond_order_text = ' '.join(bond_order_lines)

    # Now parse the bond order text
    import re
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
    for i, line in enumerate(lines):
        if 'MBIS ANALYSIS' in line:
            mbis_section = True
            continue
        if mbis_section and 'ATOM     CHARGE    POPULATION     SPIN' in line:
            # Now, read the data
            for j in range(i+1, len(lines)):
                line_content = lines[j].strip()
                if line_content == '' or 'MBIS VALENCE-SHELL DATA' in line_content or 'Total charge' in line_content:
                    break
                parts = line_content.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    atom_index = int(parts[0])
                    element = parts[1]
                    charge = float(parts[2])
                    mbis_data[atom_index] = {'element': element, 'charge': charge}
            break

    return mbis_data

def get_valence_electrons(element):
    return [a.number for a in Atoms(element)][0] - sum(1 for s in chemical_symbols[:chemical_symbols.index(element)+1] if s in ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'])

def find_electron_donors(mayer_data, mbis_data):
    donor_atoms = {}

    for atom_index, mbis_info in mbis_data.items():
        element = mbis_info['element']
        mbis_charge = mbis_info['charge']

        # Get valence electrons
        valence_electrons = get_valence_electrons(element)

        # Available electrons: valence electrons minus MBIS charge
        available_electrons = valence_electrons - mbis_charge

        # Get bonding electrons from Mayer bond orders
        bonding_electrons = sum(order for (a1, a2), order in mayer_data['bond_orders'].items() if atom_index in (a1, a2))

        # Calculate free electrons (available for external bonding)
        free_electrons = available_electrons - bonding_electrons

        # Consider atoms with significant free electrons as potential donors
        if free_electrons > 0.1:
            donor_atoms[atom_index] = {
                'element': element,
                'available_electrons': free_electrons,
                'mbis_charge': mbis_charge,
                'bonding_electrons': bonding_electrons
            }

    return donor_atoms
