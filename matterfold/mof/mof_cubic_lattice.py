# Generate cubic lattice MOF

from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import write
from ase.build import make_supercell


def cubic_lattice_mof(num_bismuth: int, ligand: str, lattice_constant=10.0, repetitions=(3, 3, 3)):
    metal_unit = generate_bismuth_polyhedron(num_bismuth)
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
    
    metal_name = f"Bi{num_bismuth}"
    ligand_name = "".join(ligand_unit.get_chemical_symbols())
    repetitions_str = f"{repetitions[0]}x{repetitions[1]}x{repetitions[2]}"
    filename = f"{metal_name}_{ligand_name}_MOF_{repetitions_str}.xyz"
    write(filename, mof)
    
    return mof
