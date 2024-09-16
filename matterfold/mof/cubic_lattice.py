# Generate cubic lattice MOF

from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import write
from ase.build import make_supercell
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones


# SMILES to RDKit to ASE Atoms object
def generate_ligand(smiles: str) -> Atoms:
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
    
    return atoms_obj

# Octohedral Bismuth Core ASE Atoms object
def generate_octahedral_bismuth() -> Atoms:
    positions = [
        (3, 0, 0),
        (-3, 0, 0),
        (0, 3, 0),
        (0, -3, 0),
        (0, 0, 3),
        (0, 0, -3)
    ]
    atoms = Atoms('Bi6', positions=positions)
    
    # lj_calc = LennardJones(sigma=3.5, epsilon=0.2)
    # atoms.calc = lj_calc
    # optimizer = BFGS(atoms)
    # optimizer.run(fmax=0.01)

    filename = 'Bi6_octahedron.xyz'
    write(filename, atoms)
    
    return atoms

# Arranges two molecules in a cubic lattice using ASE
def cubic_lattice_mof(ligand: str, lattice_constant=10.0, repetitions=(3, 3, 3)):

    metal_unit = generate_octahedral_bismuth()
    ligand_unit = generate_ligand(ligand)
    
    cell = Atoms(cell=[(lattice_constant, 0, 0), (0, lattice_constant, 0), (0, 0, lattice_constant)])
    
    metal_unit.translate([0, 0, 0])
    cell += metal_unit
    
    ligand_unit.translate([lattice_constant/2, lattice_constant/2, lattice_constant/2])
    cell += ligand_unit
    
    P = [[repetitions[0], 0, 0],
         [0, repetitions[1], 0],
         [0, 0, repetitions[2]]]
    
    mof = make_supercell(cell, P)
    
    return mof
