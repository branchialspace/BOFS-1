# utils

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Fragments
from ase import Atoms
from ase.io import write
from typing import Tuple, Dict, List, Union
import inspect
import numpy as np


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



from pymatgen.core.periodic_table import Element

# Get valence using pymatgen
def get_valence_electrons(element):
    el = Element(element)
    valence_electrons = el.full_electronic_structure
    valence = 0
    max_n = max([n for (n, l, occ) in valence_electrons])
    for (n, l, occ) in valence_electrons:
        if n == max_n or (el.is_transition_metal and n == max_n - 1):
            valence += occ
    return valence
