import os
import re
from pathlib import Path
from datetime import datetime
from pprint import pformat
import numpy as np
import spglib
from ase import Atoms
from ase.io import read, write
from ase.data import chemical_symbols
from ase.neighborlist import NeighborList, natural_cutoffs
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.structure_matcher import StructureMatcher
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
    Run QE relax/vc-relax, parse relaxed structure from output, save as new CIF.
    Returns
    relaxed_path : Path
        Path to relaxed structure CIF
    """
    # Run QE relax
    bofs1.pwx(structure_path, relax_config)
    calculation = relax_config['control']['calculation']
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    pwo_path = f"{structure_name}_{calculation}.pwo"
    with open(pwo_path, 'r') as f:
        content = f.read()
    final_block_match = re.search(
        r'Begin final coordinates.*?End final coordinates',
        content,
        re.DOTALL)
    final_block = final_block_match.group(0)
    # Parse CELL_PARAMETERS (present in vc-relax)
    cell_match = re.search(
        r'CELL_PARAMETERS \(angstrom\)\s*\n\s*([\d\.\-\+eE ]+)\n\s*([\d\.\-\+eE ]+)\n\s*([\d\.\-\+eE ]+)',
        final_block)
    if cell_match:
        cell = []
        for i in range(1, 4):
            vec = [float(x) for x in cell_match.group(i).split()]
            cell.append(vec)
    else:
        # If no new cell parameters (regular relax, not vc-relax), use original
        original_atoms = read(structure_path)
        cell = original_atoms.get_cell().tolist()
    # Parse ATOMIC_POSITIONS
    positions_match = re.search(
        r'ATOMIC_POSITIONS \(angstrom\)\s*\n(.*?)(?=End final coordinates)',
        final_block,
        re.DOTALL)
    positions_block = positions_match.group(1).strip()
    symbols = []
    positions = []
    for line in positions_block.split('\n'):
        parts = line.split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    # Create new Atoms object with relaxed coordinates
    relaxed_atoms = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=True)
    # Write new filename with abbreviated calculation right after serialization tag
    path = Path(structure_path)
    name = path.stem
    prefix = "vcr-" if calculation == "vc-relax" else "r-"
    if len(name) > 13 and name[12] == '-' and name[:12].isdigit():
        serial_tag = name[:13]
        rest = name[13:]
        new_name = f"{serial_tag}{prefix}{rest}"
    else:
        new_name = f"{prefix}{name}"
    relaxed_path = Path.cwd() / f"{new_name}.cif"
    write(relaxed_path, relaxed_atoms)
    print(f"Relaxed structure saved to: {relaxed_path}")

    return relaxed_path

def symmetrize_structure(structure_path, symprec=1e-5):
    """
    Determine space group with pymatgen.
    Write space group dataset to file with _spacegroup suffix.
    Symmetrize to primitive cell compatible with QE pw.x ibrav=0.
    """
    pmg_structure = Structure.from_file(structure_path)
    # Analyze symmetry
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
    spacegroup = sga.get_space_group_symbol()
    number = sga.get_space_group_number()
    structure_stem = os.path.splitext(structure_path)[0]
    print(f"Detected space group of structure {structure_stem}: {spacegroup} ({number})")
    # Write symmetry info
    with open(f"{structure_stem}_spacegroup", 'w') as f:
        f.write(f"Space group: {spacegroup} ({number})\n")
        f.write(f"Point group: {sga.get_point_group_symbol()}\n")
        f.write(f"Crystal system: {sga.get_crystal_system()}\n")
    # Find true primitive
    primitive_structure = sga.find_primitive()
    # Write filename
    path = Path(structure_path)
    name = path.stem
    prefix = "sym-"
    if len(name) > 13 and name[12] == '-' and name[:12].isdigit():
        serial_tag = name[:13]
        rest = name[13:]
        new_name = f"{serial_tag}{prefix}{rest}"
    else:
        new_name = f"{prefix}{name}"
    symmetrized_path = Path.cwd() / f"{new_name}.cif"
    # Write cif with pymatgen
    writer = CifWriter(primitive_structure, symprec=None, refine_struct=False)
    writer.write_file(str(symmetrized_path))
    print(f"Symmetrized structure saved to: {symmetrized_path}")

    return symmetrized_path
    
def compare_structure(structure_paths, symprec=1e-5, metal="Bi"):
    """
    Quantify geometric differences between near-identical relaxed structures.
    Compare lattice parameters, volume, and atomic displacements pairwise.
    Compare specified metal center species neighborhood topology for each structure.
    Match corresponding atoms by species using pymatgen StructureMatcher.
    Write results to file with _compare suffix.
    structure_paths : list
        List of paths to .cif files (minimum 2).
    symprec : float
        Symmetry precision for space group detection.
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
    # Analyze symmetry for each structure
    spacegroups = []
    for path, name in zip(structure_paths, names):
        pmg_structure = Structure.from_file(path)
        sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
        spacegroup = sga.get_space_group_symbol()
        number = sga.get_space_group_number()
        spacegroups.append({'symbol': spacegroup, 'number': number})
    # Compare specified metal center species neighborhood topology for each structure.
    metal_distance_stats = []
    for s, name in zip(structures, names):
        symbols = np.array(s.get_chemical_symbols())
        metal_indices = np.where(symbols == metal)[0]
        if len(metal_indices) >= 2:
            all_dist = s.get_all_distances(mic=True)
            metal_dist = all_dist[np.ix_(metal_indices, metal_indices)]
            np.fill_diagonal(metal_dist, np.inf)
            shortest = metal_dist.min()
            longest = metal_dist[metal_dist < np.inf].max()
            per_atom_avg = np.mean(metal_dist[metal_dist < np.inf].reshape(len(metal_indices), -1), axis=1)
            metal_stats = {'shortest': shortest, 'longest': longest, 'per_atom_avg': dict(zip(metal_indices, per_atom_avg))}
            print(f"\n{name} {metal}-{metal} distances:")
            print(f"  Shortest: {shortest:.4f} Å | Longest: {longest:.4f} Å")
            for idx, avg in zip(metal_indices, per_atom_avg):
                print(f"  {metal}[{idx}] avg: {avg:.4f} Å")
        else:
            metal_stats = None
            print(f"\n{name}: < 2 {metal} atoms, no {metal}-{metal} distances")
        metal_distance_stats.append(metal_stats)
    def match_atoms(s1, s2):
        """
        Match atoms between two structures using pymatgen StructureMatcher.
        Returns reordered indices for s2 to match s1 atom ordering.
        """
        symbols1 = s1.get_chemical_symbols()
        symbols2 = s2.get_chemical_symbols()
        # Validate same species composition
        if sorted(symbols1) != sorted(symbols2):
            raise ValueError("Structures have different chemical compositions")
        # Convert ASE Atoms to pymatgen Structure
        adaptor = AseAtomsAdaptor()
        pmg_s1 = adaptor.get_structure(s1)
        pmg_s2 = adaptor.get_structure(s2)
        # Use StructureMatcher to get site mapping
        # primitive_cell=False ensures we don't reduce to primitive cell
        # scale=False prevents volume rescaling
        matcher = StructureMatcher(primitive_cell=False, scale=False)
        mapping = matcher.get_mapping(pmg_s1, pmg_s2)
        if mapping is None:
            raise ValueError("StructureMatcher could not find a mapping between structures")
        # mapping[i] gives the index in s2 that corresponds to site i in s1
        reorder = np.array(mapping, dtype=int)
        return reorder
    results = {
        'names': names,
        'spacegroups': spacegroups,
        'metal_distance_stats': metal_distance_stats,
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
    compare_path = f"{names[-1]}_compare"
    with open(compare_path, 'w') as f:
        f.write(pformat(results, width=100))
        # Write metal-metal distance statistics section
        f.write(f"\n\n{'=' * 60}\n")
        f.write(f"{metal}-{metal} DISTANCE ANALYSIS\n")
        f.write(f"{'=' * 60}\n")
        for name, metal_stats in zip(names, metal_distance_stats):
            f.write(f"\n{name}:\n")
            if metal_stats:
                f.write(f"  Shortest {metal}-{metal}: {metal_stats['shortest']:.4f} Å\n")
                f.write(f"  Longest {metal}-{metal}:  {metal_stats['longest']:.4f} Å\n")
                for idx, avg in metal_stats['per_atom_avg'].items():
                    f.write(f"  {metal}[{idx}] avg: {avg:.4f} Å\n")
            else:
                f.write(f"  < 2 {metal} atoms\n")
        # Write full CIF files for reference
        for i, (name, path) in enumerate(zip(names, structure_paths), start=1):
            abs_path = Path(path).absolute()
            f.write(f"\n\n{'=' * 60}\n")
            f.write(f"[{i}] {name}\n")
            f.write(f"    {abs_path}\n")
            f.write(f"{'=' * 60}\n")
            with open(path, 'r') as cif:
                f.write(cif.read())
    print(f"\nComparison results written to: {compare_path}")

    return results
