import os
from pathlib import Path
from datetime import datetime
import spglib
from ase import Atoms
from ase.io import read, write
import bofs1
from optimade.client import OptimadeClient
from optimade.adapters import Structure as OptimadeStructure


def get_structure(source, material_id):
    """
    Fetches a structure from an OPTIMADE provider, saves it as a CIF, 
    and returns the absolute file path.
    source : str
    	The database alias ('mp', 'cod', 'nomad', 'oqmd', 'aflow') or a full valid OPTIMADE URL.
    material_id : str
    	The ID of the material (e.g., 'mp-23152' or '1010068').
    Returns 
    structure_path : str
    	Path to the saved structure.
    """
    provider_urls = {
        "mp": "https://optimade.materialsproject.org",
        "cod": "https://www.crystallography.net/cod/optimade",
        "oqmd": "http://oqmd.org/optimade/",
        "nomad": "https://nomad-lab.eu/prod/rae/optimade/",
        "aflow": "http://aflow.org/API/optimade/"}
    base_url = provider_urls.get(source, source)
    client = OptimadeClient(base_urls=[base_url])
    print(f"Querying {base_url} for {material_id}...")
    filter_str = f'id="{material_id}"'
    client.get(filter=filter_str)
    results = client.all_results["structures"][filter_str][base_url]["data"]
    entry = results[0]
    pmg_structure = OptimadeStructure(entry).as_pymatgen
    structure_path = os.path.abspath(f"{material_id}.cif")
    pmg_structure.to(filename=structure_path)
    print(f"Structure saved to: {structure_path}")
    
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


