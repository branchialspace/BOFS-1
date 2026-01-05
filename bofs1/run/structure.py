import os
from pathlib import Path
from datetime import datetime
import numpy as np
from scipy.optimize import linear_sum_assignment
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

def compare_structure(structure_paths):
    """
    Quantify geometric differences between near-identical relaxed structures.
    Compares lattice parameters, volume, and atomic displacements pairwise.
    Matches corresponding atoms by species using Hungarian algorithm.
    Write results to file with _compare suffix.
    structure_paths : list
        List of paths to .cif files (minimum 2).
    Returns
    results : dict
        Dictionary containing pairwise comparison metrics.
    """
    if len(structure_paths) < 2:
        raise ValueError("At least 2 structure paths required for comparison.")
    structures = [read(p) for p in structure_paths]
    names = [Path(p).stem for p in structure_paths]
    n_structures = len(structures)
    # Validate atom counts match
    n_atoms = len(structures[0])
    for i, s in enumerate(structures):
        if len(s) != n_atoms:
            raise ValueError(f"Atom count mismatch: {names[0]} has {n_atoms}, {names[i]} has {len(s)}")
    def match_atoms(s1, s2):
        """
        Match atoms between two structures by species using Hungarian.
        Minimizes total displacement with periodic boundary conditions.
        Returns reordered indices for s2 to match s1 atom ordering.
        """
        symbols1 = s1.get_chemical_symbols()
        symbols2 = s2.get_chemical_symbols()
        frac1 = s1.get_scaled_positions()
        frac2 = s2.get_scaled_positions()
        # Validate same species composition
        if sorted(symbols1) != sorted(symbols2):
            raise ValueError("Structures have different chemical compositions")
        reorder = np.zeros(len(s2), dtype=int)
        species = sorted(set(symbols1))
        for sp in species:
            # Indices of this species in each structure
            idx1 = [i for i, sym in enumerate(symbols1) if sym == sp]
            idx2 = [i for i, sym in enumerate(symbols2) if sym == sp]
            if len(idx1) != len(idx2):
                raise ValueError(f"Species count mismatch for {sp}")
            # Fractional positions for this species
            pos1 = frac1[idx1]
            pos2 = frac2[idx2]
            # Cost matrix: minimum image distance between all pairs
            n_sp = len(idx1)
            cost = np.zeros((n_sp, n_sp))
            for a in range(n_sp):
                diff = pos2 - pos1[a]
                diff = diff - np.round(diff)
                cost[a, :] = np.linalg.norm(diff, axis=1)
            # Hungarian algorithm for optimal assignment
            row_ind, col_ind = linear_sum_assignment(cost)
            # Map s1 indices to matched s2 indices
            for r, c in zip(row_ind, col_ind):
                reorder[idx1[r]] = idx2[c]
        return reorder
    results = {
        'names': names,
        'pairwise': []}
    for i in range(n_structures):
        for j in range(i + 1, n_structures):
            s1, s2 = structures[i], structures[j]
            name1, name2 = names[i], names[j]
            # Match corresponding atoms
            reorder = match_atoms(s1, s2)
            # Lattice parameters [a, b, c, alpha, beta, gamma]
            cell1 = np.array(s1.cell.cellpar())
            cell2 = np.array(s2.cell.cellpar())
            cell_diff = cell2 - cell1
            cell_pct = 100.0 * cell_diff / np.where(cell1 != 0, cell1, 1.0)
            # Volume
            vol1, vol2 = s1.get_volume(), s2.get_volume()
            vol_diff = vol2 - vol1
            vol_pct = 100.0 * vol_diff / vol1 if vol1 != 0 else 0.0
            # Fractional position displacements with matched atoms
            frac1 = s1.get_scaled_positions()
            frac2 = s2.get_scaled_positions()[reorder]
            frac_disp = frac2 - frac1
            # Wrap to [-0.5, 0.5] for periodic boundary
            frac_disp = frac_disp - np.round(frac_disp)
            # Convert to Cartesian using average cell
            avg_cell = 0.5 * (np.array(s1.get_cell()) + np.array(s2.get_cell()))
            cart_disp = frac_disp @ avg_cell
            # Displacement magnitudes
            disp_mag = np.linalg.norm(cart_disp, axis=1)
            rmsd = np.sqrt(np.mean(disp_mag ** 2))
            max_disp = np.max(disp_mag)
            mean_disp = np.mean(disp_mag)
            max_disp_idx = int(np.argmax(disp_mag))
            max_disp_atom = s1.get_chemical_symbols()[max_disp_idx]
            # Per-species RMSD
            symbols = s1.get_chemical_symbols()
            species = sorted(set(symbols))
            species_rmsd = {}
            for sp in species:
                mask = np.array([sym == sp for sym in symbols])
                if np.any(mask):
                    species_rmsd[sp] = np.sqrt(np.mean(disp_mag[mask] ** 2))
            comparison = {
                'pair': (name1, name2),
                'cell_params_1': cell1.tolist(),
                'cell_params_2': cell2.tolist(),
                'cell_diff': cell_diff.tolist(),
                'cell_pct_change': cell_pct.tolist(),
                'volume_1': vol1,
                'volume_2': vol2,
                'volume_diff': vol_diff,
                'volume_pct_change': vol_pct,
                'rmsd': rmsd,
                'mean_displacement': mean_disp,
                'max_displacement': max_disp,
                'max_displacement_atom': (max_disp_idx, max_disp_atom),
                'species_rmsd': species_rmsd,
                'displacements': disp_mag.tolist()}
            results['pairwise'].append(comparison)
            print(f"\n{name1} → {name2}:")
            print(f"  Volume: {vol1:.3f} → {vol2:.3f} Å³ ({vol_pct:+.2f}%)")
            print(f"  RMSD: {rmsd:.4f} Å | Max: {max_disp:.4f} Å ({max_disp_atom}[{max_disp_idx}])")
            print(f"  Cell: Δa={cell_diff[0]:.4f} Δb={cell_diff[1]:.4f} Δc={cell_diff[2]:.4f} Å")
    # Write results to file
    compare_path = f"{names[0]}_compare"
    with open(compare_path, 'w') as f:
        f.write(str(results))
    print(f"\nComparison results written to: {compare_path}")
    
    return results
