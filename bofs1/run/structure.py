import os
from pathlib import Path
from datetime import datetime
import spglib
from ase import Atoms
from ase.io import read, write
import bofs1
from pymatgen.ext.optimade import OptimadeRester

def get_structure(source, material_id):
    """
    Download a structure from a specific OPTIMADE provider, save it 
    to the current working directory, return the path.
    source : str
    	The database alias (ex: 'mp', 'cod', 'oqmd', 'aflow').
    material_id : str
    	The ID of the material (e.g., 'mp-149' or '1010068').
    Returns 
    structure_path : str
    	Path to the saved structure.
    """
    with OptimadeRester(source) as optimade:
        print(f"Querying {source.upper()} for {material_id}...")
        results = optimade.get_structures(filter=f'id="{material_id}"')
        structure_path = os.path.abspath(f"{material_id}.cif")
        results[material_id].to(filename=structure_path)
        print(f"{material_id} saved as: {structure_path}")
        
    return structure_path

def normalize_structure(structure_path, relax_config):
    """
    Add timestamp to structure filename.
    Run QE vc-relax, update structure.
    Determine spacegroup with spglib.
    Returns
    relaxed_structure_path : str
        Path to relaxed structure
    """
    atoms = read(structure_path)
    # Serialize
    path = Path(structure_path)
    serial_name = f"{datetime.now():%Y%m%d%H%M}-{path.name}"
    serial_path = Path.cwd() / serial_name
    write(serial_path, atoms)
    # Run QE vc-relax
    bofs1.pwx(serial_path, relax_config)
    relaxed_atoms = read(serial_path)
    # Detect symmetry with spglib
    lattice = relaxed_atoms.get_cell()
    positions = relaxed_atoms.get_scaled_positions()
    numbers = relaxed_atoms.get_atomic_numbers()
    cell_tuple = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell_tuple)
    if isinstance(dataset, dict):
        spacegroup = dataset["international"]
        number = dataset["number"]
    else:
        spacegroup = dataset.international
        number = dataset.number
    print(f"Detected space group after vc-relax: {spacegroup} ({number})")

    return str(serial_path)


