# Mayer and Lowdin electron-donor analysis in ORCA

import os
from ase import Atoms
from ase.calculators.orca import OrcaProfile, ORCA
from ase.data import chemical_symbols

def electron_donors(ligand, charge=0, mult=1, method='B3LYP', basis_set='def2-SVP'):
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
    Print[P_Mayer] 1
    Print[P_Loewdin] 1
    Print[P_AtCharges_L] 1
    Print[P_BondOrder_L] 1
    Print[P_AtPopMO_L] 1
    Print[P_OrbPopMO_L] 1
end

%method
    MAYER_BONDORDERTHRESH 0.05
    LOEWDIN_BONDORDERTHRESH 0.05
end
"""
    )
    
    ligand.calc = calculator
    try:
        energy = ligand.get_potential_energy()
        print(f"Calculation completed. Energy: {energy}")
    except Exception as e:
        print(f"Error during ORCA calculation: {e}")
    
    orca_output_file = 'orca.out'
    #mayer_data = parse_mayer_data(orca_output_file)
    #lowdin_data = parse_lowdin_data(orca_output_file)
    #donor_atoms = find_electron_donors(mayer_data, lowdin_data)
    
    return orca_output_file #donor_atoms

def parse_mayer_data(filename):
    mayer_data = {'charges': {}, 'bond_orders': {}, 'free_valences': {}}
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse Mayer charges and free valences
    mayer_section = False
    for i, line in enumerate(lines):
        if 'MAYER POPULATION ANALYSIS' in line:
            mayer_section = True
            continue
        if mayer_section and 'ATOM' in line and 'NA' in line:
            for j in range(i+2, len(lines)):
                parts = lines[j].strip().split()
                if len(parts) < 7:
                    break
                atom_index = int(parts[0])
                element = parts[1]
                charge = float(parts[3])  # QA - Mulliken gross atomic charge
                free_valence = float(parts[6])  # FA - Mayer's free valence
                mayer_data['charges'][atom_index] = {'element': element, 'charge': charge}
                mayer_data['free_valences'][atom_index] = free_valence
            break

    # Parse Mayer bond orders
    bond_order_section = False
    for line in lines:
        if 'Mayer bond orders larger than' in line:
            bond_order_section = True
            continue
        if bond_order_section:
            if line.strip() == '':
                break
            parts = line.strip().split()
            i = 0
            while i < len(parts):
                if parts[i].startswith('B('):
                    atom_pair = parts[i].split('(')[1].split(')')[0]
                    atom1, atom2 = map(lambda x: int(x.strip('-')), atom_pair.split(','))
                    order = float(parts[i+2])
                    mayer_data['bond_orders'][(atom1, atom2)] = order
                    i += 3
                else:
                    i += 1

    return mayer_data

def parse_lowdin_data(filename):
    lowdin_data = {'charges': {}, 'bond_orders': {}, 'orbital_populations': {}}
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse Löwdin charges
    charge_section = False
    for line in lines:
        if 'LOEWDIN ATOMIC CHARGES' in line:
            charge_section = True
            continue
        if charge_section:
            if line.strip() == '':
                charge_section = False
                break
            parts = line.strip().split()
            if len(parts) == 3:
                atom_index = int(parts[0])
                element = parts[1].rstrip(':')
                charge = float(parts[2])
                lowdin_data['charges'][atom_index] = {'element': element, 'charge': charge}

    # Parse Löwdin bond orders
    bond_order_section = False
    for line in lines:
        if 'LOEWDIN BOND ORDERS (THRESH' in line:
            bond_order_section = True
            continue
        if bond_order_section:
            if line.strip() == '':
                bond_order_section = False
                break
            parts = line.strip().split()
            for i in range(0, len(parts), 6):
                if i + 5 < len(parts):
                    atom1 = int(parts[i+1].strip('()-,'))
                    atom2 = int(parts[i+3].strip('()-,'))
                    order = float(parts[i+5].strip('()-,'))
                    lowdin_data['bond_orders'][(atom1, atom2)] = order

    # Parse Löwdin orbital populations
    orbital_pop_section = False
    current_atom = None
    for i, line in enumerate(lines):
        if 'LOEWDIN ORBITAL POPULATIONS PER MO' in line:
            orbital_pop_section = True
            continue
        if orbital_pop_section:
            if 'LOEWDIN REDUCED ORBITAL CHARGES' in line:
                break
            if line.startswith('THRESHOLD FOR PRINTING IS'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1] in chemical_symbols:
                current_atom = int(parts[0])
                lowdin_data['orbital_populations'][current_atom] = {}
            elif current_atom is not None and len(parts) >= 2:
                orbital = parts[0] + parts[1]
                populations = [float(p) for p in parts[2:] if p != '---------']
                lowdin_data['orbital_populations'][current_atom][orbital] = sum(populations)

    return lowdin_data

def find_electron_donors(mayer_data, lowdin_data):
    donor_atoms = {}
    
    for atom_index, mayer_info in mayer_data['charges'].items():
        element = mayer_info['element']
        mayer_charge = mayer_info['charge']
        free_valence = mayer_data['free_valences'][atom_index]

        # Use Mayer free valence as the primary indicator of available electrons
        if free_valence > 0.1:
            # Calculate bonding electrons using Mayer bond orders
            bonding_electrons = sum(order for (a1, a2), order in mayer_data['bond_orders'].items() if atom_index in (a1, a2))
            
            # Use Löwdin orbital populations for a more detailed electron count
            valence_electrons = get_valence_electrons(element)
            if atom_index in lowdin_data['orbital_populations']:
                available_electrons = sum(pop for orbital, pop in lowdin_data['orbital_populations'][atom_index].items() if orbital.lower() not in ['1s', '2s', '2p'])
            else:
                available_electrons = valence_electrons - bonding_electrons + mayer_charge

            donor_atoms[atom_index] = {
                'element': element,
                'available_electrons': available_electrons,
                'free_valence': free_valence,
                'mayer_charge': mayer_charge
            }

    return donor_atoms

def get_valence_electrons(element):
    return [a.number for a in Atoms(element)][0] - sum(1 for s in chemical_symbols[:chemical_symbols.index(element)+1] if s in ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'])
