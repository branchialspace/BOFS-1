import os
import re
from pathlib import Path
from datetime import datetime
import spglib
from ase import Atoms
from ase.io import read, write

def normalize_structure(structure_path):
    """
    Determine spacegroup of structure and write symmetrized unit cell.
    Add timestamp to structure filename.
    Returns
    norm_structure : symmetrized cif 
    """
    atoms = read(structure_path)
    # Symmetrize (spglib)
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell_tuple = (lattice, positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell_tuple)
    spacegroup = dataset['international'] if isinstance(dataset, dict) else dataset.international
    number = dataset['number'] if isinstance(dataset, dict) else dataset.number
    print(f"Detected space group: {spacegroup} ({number})")
    standardized_cell = spglib.standardize_cell(
        cell_tuple,
        to_primitive=False,
        no_idealize=False)
    std_lattice, std_positions, std_numbers = standardized_cell
    atoms_std = Atoms(
        numbers=std_numbers,
        cell=std_lattice,
        scaled_positions=std_positions,
        pbc=True)
    # Serialize
    path = Path(structure_path)
    serial_name = f"{datetime.now():%Y%m%d%H%M}-{path.name}"
    serial_path = Path.cwd() / serial_name
    write(serial_path, atoms_std)
    
    return str(serial_path)
