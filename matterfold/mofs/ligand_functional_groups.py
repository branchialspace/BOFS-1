# Identify the functional groups of a given ligand

import numpy as np
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Fragments
from rdkit.Chem import AllChem


def ase_to_rdkit(atoms):
    # Get atomic numbers and coordinates from ASE Atoms object
    atomic_numbers = atoms.get_atomic_numbers()
    positions = atoms.get_positions()
    
    # Create an RDKit molecule
    mol = Chem.Mol()
    editable_mol = Chem.EditableMol(mol)
    
    # Add atoms to the molecule
    for atomic_num, pos in zip(atomic_numbers, positions):
        atom = Chem.Atom(int(atomic_num))
        atom_idx = editable_mol.AddAtom(atom)
        conf = editable_mol.GetMol().GetConformer()
        conf.SetAtomPosition(atom_idx, pos)
    
    # Convert editable molecule to a regular molecule
    mol = editable_mol.GetMol()
    
    # Sanitize the molecule
    Chem.SanitizeMol(mol)
    
    return mol

def identify_functional_groups(atoms):
    # Convert ASE Atoms to RDKit molecule
    mol = ase_to_rdkit(atoms)
    
    functional_groups = []
    
    # Use rdMolDescriptors to get functional group counts
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    fg_counts = rdMolDescriptors.GetMorganFeatureFingerprint(mol, useFeatures=True, nBits=2048)
    
    # Iterate through all available Fragment functions
    for name in dir(Fragments):
        if name.startswith('fr_'):
            func = getattr(Fragments, name)
            count = func(mol)
            if count > 0:
                functional_groups.append((name[3:], count))
    
    return functional_groups
