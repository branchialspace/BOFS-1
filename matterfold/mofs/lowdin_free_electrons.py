# Lowdin free electron analysis in ORCA

import os
from ase import Atoms
from ase.calculators.orca import OrcaProfile, ORCA
from ase.data import chemical_symbols


def lowdin_electron_donors(ligand, charge=0, mult=1, method='B3LYP', basis_set='def2-SVP'):
    # Define the ORCA calculator
    orca_path = '/root/orca_6_0_0/orca'
    profile = OrcaProfile(command=orca_path)

    calculator = ORCA(
        profile=profile,
        label='orca_calculation',
        task='SP',
        charge=charge,
        mult=mult,
        orcasimpleinput=f'{method} {basis_set}',
        orcablocks="""
%output
    Print[P_Loewdin] 1
    Print[P_AtCharges_L] 1
    Print[P_BondOrder_L] 1
end
"""
    )
    
    ligand.calc = calculator
    try:
        energy = ligand.get_potential_energy()
        print(f"Calculation completed. Energy: {energy}")
    except Exception as e:
        print(f"Error during ORCA calculation: {e}")
    
    # Parse the ORCA output file
    orca_output_file = 'orca.out'
    
    lowdin_data = parse_lowdin_data(orca_output_file)
    
    # Analyze results
    donor_atoms = find_free_electrons(lowdin_data)
    
    return donor_atoms

# Parse the ORCA output file to extract Löwdin charges and bond orders
def parse_lowdin_data(filename):
    lowdin_data = {'charges': {}, 'bond_orders': {}}
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse Löwdin charges
    charge_section = False
    for line in lines:
        if 'LOEWDIN ATOMIC CHARGES' in line:
            charge_section = True
            continue
        if charge_section:
            if '----------' in line:
                charge_section = False
                break
            parts = line.strip().split()
            if len(parts) >= 3:
                atom_index = int(parts[0])
                element = parts[1]
                charge = float(parts[2])
                lowdin_data['charges'][atom_index] = {'element': element, 'charge': charge}

    # Parse Löwdin bond orders
    bond_order_section = False
    for line in lines:
        if 'Loewdin Bond Orders' in line:
            bond_order_section = True
            continue
        if bond_order_section:
            if '--------' in line:
                break
            parts = line.strip().split()
            if len(parts) >= 3:
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                order = float(parts[2])
                lowdin_data['bond_orders'][(atom1, atom2)] = order

    return lowdin_data

# Analyze the electron distribution to identify electron-donor atoms and their available free electrons
def find_free_electrons(lowdin_data):
    donor_atoms = {}
    lowdin_charges = lowdin_data['charges']
    bond_orders = lowdin_data['bond_orders']

    for atom_index, lowdin_info in lowdin_charges.items():
        element = lowdin_info['element']
        lowdin_charge = lowdin_info['charge']

        # Identify potential donors based on negative charge
        if lowdin_charge < -0.1:
            # Calculate free electrons
            valence_electrons = get_valence_electrons(element)
            bonding_electrons = sum(order for (a1, a2), order in bond_orders.items() if atom_index in (a1, a2))
            
            # Free electrons = valence electrons - bonding electrons + excess negative charge
            free_electrons = valence_electrons - bonding_electrons - lowdin_charge

            if free_electrons > 0:
                donor_atoms[atom_index] = {
                    'element': element,
                    'free_electrons': free_electrons,
                    'lowdin_charge': lowdin_charge
                }

    return donor_atoms

# Return the number of valence electrons for an element using ASE
def get_valence_electrons(element):
    return [a.number for a in Atoms(element)][0] - sum(1 for s in chemical_symbols[:chemical_symbols.index(element)+1] if s in ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'])
