# Generate ligand, identify functional groups

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Fragments
from ase import Atoms
from ase.io import write
from typing import Tuple, Dict, List, Union
import inspect
import numpy as np


# Generate ligand: SMILES to RDKit to ASE Atoms object
def generate_ligand(smiles: str) -> Tuple[Atoms, Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    
    conformer = mol.GetConformer()
    atoms = []
    positions = []
    
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atoms.append(symbol)
    
    for i in range(mol.GetNumAtoms()):
        pos = conformer.GetAtomPosition(i)
        positions.append([pos.x, pos.y, pos.z])
    
    atoms_obj = Atoms(atoms, positions=positions)
    filename = "".join(atoms_obj.get_chemical_symbols()) + ".xyz"
    write(filename, atoms_obj)
    
    return atoms_obj, mol

# Identify functional groups from RDKit Fragments    
def identify_functional_groups(mol: Chem.Mol) -> Dict[int, Dict[str, Union[str, List[int]]]]:
    # Sanitize the molecule to address aromatic atom issues
    Chem.SanitizeMol(mol)
    
    functional_groups = {}
    instance_counter = 0
    
    # Get all fragment functions from rdkit.Chem.Fragments
    fragment_functions = [f for f in dir(Fragments) if f.startswith('fr_')]
    
    for func_name in fragment_functions:
        func = getattr(Fragments, func_name)
        sig = inspect.signature(func)
        pattern = sig.parameters['pattern'].default
        if pattern is None:
            continue  # Skip functions without a default 'pattern' parameter
        try:
            smarts_mol = Chem.MolFromSmarts(pattern)
            if smarts_mol is None:
                continue  # Invalid SMARTS pattern
            matches = mol.GetSubstructMatches(smarts_mol)
        except:
            continue  # In case of any issues with the SMARTS pattern
        
        # Get the SMARTS pattern for this functional group
        fg_smarts = Chem.MolToSmarts(smarts_mol)
                        
        # For each match, create a new entry in the dictionary
        for match in matches:
            functional_groups[instance_counter] = {
                "type": func_name,
                "smarts": fg_smarts,
                "atom_indices": list(match)
            }
            instance_counter += 1
    
    return functional_groups

# Calculate partial charges using Gasteiger method
def calculate_partial_charges(mol: Chem.Mol) -> List[float]:
    AllChem.ComputeGasteigerCharges(mol)
    charges = []
    for atom in mol.GetAtoms():
        try:
            charge = float(atom.GetProp('_GasteigerCharge'))
        except:
            charge = 0.0  # Default to 0 if charge not available
        charges.append(charge)
    return charges

# Identify donor atoms (N, O, P, S with lone pairs and non-positive formal charge)
def identify_donor_atoms(mol: Chem.Mol) -> List[int]:
    donor_atoms = []
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num in [7, 8, 15, 16]:  # N, O, P, S
            if atom.GetFormalCharge() <= 0:
                donor_atoms.append(atom.GetIdx())
    return donor_atoms

# Calculate available electrons for donation
def calculate_available_electrons(mol: Chem.Mol, donor_atoms: List[int], charges: List[float]) -> Dict[str, Union[float, Dict[int, float]]]:
    # Define base lone pair electrons for each donor atom type
    base_lone_pairs = {
        7: 2,   # Nitrogen
        8: 2,   # Oxygen
        15: 2,  # Phosphorus
        16: 2   # Sulfur
    }
    
    available_electrons_per_atom = {}
    total_available_electrons = 0.0
    
    for idx in donor_atoms:
        atom = mol.GetAtomWithIdx(idx)
        atomic_num = atom.GetAtomicNum()
        symbol = atom.GetSymbol()
        base_electrons = base_lone_pairs.get(atomic_num, 0)
        partial_charge = charges[idx]
        
        # Charge adjustment: more negative charge increases available electrons
        # Here, we assume that each -1 charge adds 1 electron, scaled by partial charge
        charge_adjustment = -partial_charge  # Negative partial charge increases electrons
        
        # Available electrons calculation
        available_electrons = base_electrons + charge_adjustment
        
        # Ensure available electrons are within reasonable bounds
        available_electrons = max(base_electrons, available_electrons)  # At least base electrons
        available_electrons = min(available_electrons, base_electrons + 2)  # Max 2 extra electrons
        
        # Round to two decimal places for clarity
        available_electrons = round(available_electrons, 2)
        
        available_electrons_per_atom[idx] = available_electrons
        total_available_electrons += available_electrons
    
    # Round total to two decimal places
    total_available_electrons = round(total_available_electrons, 2)
    
    result = {
        'total_available_electrons': total_available_electrons,
        'available_electrons_per_atom': available_electrons_per_atom
    }
    
    return result

def analyze_ligand(smiles: str) -> Dict[str, Union[Atoms, Chem.Mol, List[float], Dict, List, float, Dict[int, float]]]:
    # Generate the ligand
    atoms_obj, mol = generate_ligand(smiles)
    
    # Identify functional groups
    functional_groups = identify_functional_groups(mol)
    
    # Calculate partial charges
    charges = calculate_partial_charges(mol)
    
    # Identify donor atoms
    donor_atoms = identify_donor_atoms(mol)
    
    # Calculate available electrons for donation
    electrons_info = calculate_available_electrons(mol, donor_atoms, charges)
    
    # Compile results
    result = {
        'atoms_obj': atoms_obj,
        'mol': mol,
        'charges': charges,
        'functional_groups': functional_groups,
        'donor_atoms': donor_atoms,
        'total_available_electrons': electrons_info['total_available_electrons'],
        'available_electrons_per_atom': electrons_info['available_electrons_per_atom']
    }

    print("Ligand Atoms with Partial Charges:")
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        charge = charges[idx]
        print(f"Atom {idx}: {symbol}, Charge: {charge:.2f}")
    
    print("\nFunctional Groups Identified:")
    for idx, fg in functional_groups.items():
        print(f"Functional Group {idx}:")
        print(f"  Type: {fg['type']}")
        print(f"  SMARTS: {fg['smarts']}")
        print(f"  Atom Indices: {fg['atom_indices']}")
    
    print("\nDonor Atoms Identified and Available Electrons:")
    for idx in donor_atoms:
        atom = mol.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        charge = charges[idx]
        available_electrons = electrons_info['available_electrons_per_atom'][idx]
        print(f"Atom {idx}: {symbol}, Charge: {charge:.2f}, Available Electrons: {available_electrons}")
    
    print(f"\nTotal Available Electrons for Donation: {electrons_info['total_available_electrons']}")
    
    return result
