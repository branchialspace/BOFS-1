# Generate cubic lattice MOF

from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.build import make_supercell


# Convert SMILES string to ASE Atoms object
def generate_mols(smiles: str) -> Atoms:
  
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
    
    return Atoms(atoms, positions=positions)

# Arranges two molecules in a cubic lattice using ASE
def cubic_lattice_mof(metal: str, ligand: str, lattice_constant=10.0, repetitions=(3, 3, 3)):

    metal_unit = generate_mols(metal)
    ligand_unit = generate_mols(ligand)
    
    cell = Atoms(cell=[(lattice_constant, 0, 0), (0, lattice_constant, 0), (0, 0, lattice_constant)])
    
    metal_unit.translate([0, 0, 0])
    cell += ase_mol1
    
    ligand_unit.translate([lattice_constant/2, lattice_constant/2, lattice_constant/2])
    cell += ase_mol2
    
    P = [[repetitions[0], 0, 0],
         [0, repetitions[1], 0],
         [0, 0, repetitions[2]]]
    
    mof = make_supercell(cell, P)
    
    return mof
