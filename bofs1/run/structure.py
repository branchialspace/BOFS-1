import os
from pathlib import Path
from datetime import datetime
import spglib
from ase import Atoms
from ase.io import read, write
import bofs1
from mp_api.client import MPRester


def get_structure(mp_api_key, material_id):
    """
    Fetches a structure from the Materials Project API, saves it as a CIF, 
    and returns the absolute file path.
    mp_api_key : str
        Materials Project API key.
    material_id : str
    	The ID of the material (e.g., 'mp-23152').
    Returns
    structure_path : str
    	Path to the saved structure.
    """
    with MPRester(mp_api_key) as mpr:
        print(f"Querying Materials Project for {material_id}...")
        pmg_structure = mpr.get_structure_by_material_id(material_id)
    structure_path = os.path.abspath(f"{material_id}.cif")
    pmg_structure.to(filename=structure_path)
    print(f"Structure saved to: {structure_path}")
    
    return structure_path

def serialize_structure(structure_path):
    """
    Add timestamp to structure filename.
    Overwrites serialization if present.
    Returns
    serial_path : Path
        Path to serialized structure
    """
    atoms = read(structure_path)
    path = Path(structure_path)
    name = path.name
    if name[12:13] == '-' and name[:12].isdigit():
        name = name[13:]
    serial_name = f"{datetime.now():%Y%m%d%H%M}-{name}"
    serial_path = Path.cwd() / serial_name
    write(serial_path, atoms)
    
    return serial_path

def relax_structure(structure_path, relax_config):
    """
    Run QE vc-relax, save relaxed structure.
    Returns
    relaxed_path : Path
        Path to relaxed structure
    """
    # Run QE vc-relax
    bofs1.pwx(structure_path, relax_config)
    
    return structure_path
    
def spglib_structure(structure_path):
    """
    Determine space group with spglib.
    Write dataset to file with _spglib suffix.
    """
    atoms = read(structure_path)
    # Detect symmetry with spglib
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell_tuple = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell_tuple)
    if isinstance(dataset, dict):
        spacegroup = dataset["international"]
        number = dataset["number"]
    else:
        spacegroup = dataset.international
        number = dataset.number
    print(f"Detected space group: {spacegroup} ({number})")
    # Write spglib dataset to file
    spglib_path = f"{serial_path}_spglib"
    with open(spglib_path, 'w') as f:
        f.write(str(dataset))
