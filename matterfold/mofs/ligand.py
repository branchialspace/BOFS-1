# Generate ligand, identify functional groups

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Fragments
from ase import Atoms
from ase.io import write
from typing import Tuple, Dict, List, Union
import inspect


# Generate ligand: SMILES to RDKit to ASE Atoms object
def generate_ligand(smiles: str) -> Tuple[Atoms, Chem.Mol]:
    mol = Chem.MolFromSmiles(smiles)
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
def identify_functional_groups(smiles: str) -> Dict[int, Dict[str, Union[str, List[int]]]]:
_, mol = generate_ligand(smiles)

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
    matches = mol.GetSubstructMatches(pattern)
    
    # Get the SMARTS pattern for this functional group
    fg_smarts = Chem.MolToSmarts(pattern)
                    
    # For each match, create a new entry in the dictionary
    for match in matches:
        functional_groups[instance_counter] = {
            "type": func_name,
            "smarts": fg_smarts,
            "atom_indices": list(match)
        }
        
        # Print the information for this instance
        print(f"Functional Group {instance_counter}:")
        print(f"  Type: {func_name}")
        print(f"  SMARTS: {fg_smarts}")
        print(f"  Atom Indices: {list(match)}")
        
        instance_counter += 1

return functional_groups
