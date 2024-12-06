# MOF cubic lattice nanoparticle

from ase import Atoms
from ase.io import write
from ase.build.supercells import find_optimal_cell_shape, make_supercell
from scipy.spatial import ConvexHull
import numpy as np


def mof_supercell(
    unit_cell: Atoms,
    metal_center: Atoms,
    ligand: Atoms,
    combined_structure: Atoms,
    bonding_sites: list,
    target_size: int = 8, # 27, 64, 125
    target_shape: str = 'sc',
    wrap: bool = False
) -> Atoms:
    """
    Creates an extended cubic lattice nanoparticle from a MOF unit cell and removes
    surface metal atoms and ligands along the convex hull.

    Parameters:
    - unit_cell: ASE Atoms object of the MOF primitive cubic unit cell
    - metal_center: ASE Atoms object of the metal cluster
    - ligand: ASE Atoms object of the ligand
    - combined_structure: ASE Atoms object of the combined metal and ligand structure
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site
    - target_size: Target number of unit cells in the supercell (default: 8 for 2x2x2)
    - target_shape: Desired supercell shape, 'sc' for simple cubic or 'fcc' for face-centered cubic (default: 'sc')
    - wrap: Whether to wrap atoms at the boundaries (default: False)

    Returns:
    - extended_lattice: ASE Atoms object of the extended MOF lattice with surface atoms and ligands removed
    """
    n_metal_center = len(metal_center)
    n_ligand = len(ligand)
    n_total = len(unit_cell)

    # Create the supercell using the transformation matrix
    P = find_optimal_cell_shape(
        cell=unit_cell.cell,
        target_size=target_size,
        target_shape=target_shape,
        verbose=True
    )

    extended_lattice = make_supercell(
        prim=unit_cell,
        P=P,
        wrap=wrap,
        order="cell-major"  # Keep atoms from same unit cell together
    )

    # Remove surface structures to maintain charge balance and geometry for aperiodic structure

    # Compute the convex hull of the extended lattice
    positions = extended_lattice.get_positions()
    hull = ConvexHull(positions)

    # Each row of hull.equations is of form [A, B, C, D] defining a plane A*x + B*y + C*z + D = 0
    equations = hull.equations
    tol = 1e-2
    on_hull_indices = []

    for i, pos in enumerate(positions):
        # Check each facet inequality:
        # For a point inside the hull, (A*x+B*y+C*z+D) < 0 for all facets (if normals point outward).
        # Points on or outside will satisfy (A*x+B*y+C*z+D) >= -tol for at least one facet.
        distances = np.dot(equations[:, :3], pos) + equations[:, 3]
        # If this point is on or outside one of the hull facets (within a small tolerance),
        # consider it a hull atom.
        if np.any(distances >= -tol):
            on_hull_indices.append(i)

    hull_indices = np.array(on_hull_indices)

    # Initialize a set to collect indices of atoms to remove
    atoms_to_remove = set()

    # Process atoms on the convex hull
    for idx in hull_indices:
        # Determine which unit cell this atom belongs to
        unit_cell_index = idx // n_total
        atom_in_unit_cell_index = idx % n_total

        if atom_in_unit_cell_index < n_metal_center:
            # It's a metal center atom
            # Remove this atom
            atoms_to_remove.add(idx)
        else:
            # It's a ligand atom
            # Remove only the ligand connected to this atom
            # First, find all ligand atoms connected to this atom in the unit cell
            ligand_start_idx = unit_cell_index * n_total + n_metal_center
            ligand_end_idx = unit_cell_index * n_total + n_total
            ligand_atom_indices = list(range(ligand_start_idx, ligand_end_idx))

            # Remove only the ligand atoms corresponding to this ligand molecule
            # Assuming ligands are added in the same order and structure as in the ligand Atoms object
            n_atoms_per_ligand = len(ligand)
            # Find which ligand molecule this atom belongs to
            ligand_molecule_index = (atom_in_unit_cell_index - n_metal_center) // n_atoms_per_ligand
            ligand_molecule_start = ligand_start_idx + ligand_molecule_index * n_atoms_per_ligand
            ligand_molecule_end = ligand_molecule_start + n_atoms_per_ligand
            ligand_molecule_indices = list(range(ligand_molecule_start, ligand_molecule_end))
            atoms_to_remove.update(ligand_molecule_indices)

            # Find the closest metal center atom to this ligand molecule
            ligand_positions = extended_lattice.get_positions()[ligand_molecule_indices]
            metal_center_indices = [unit_cell_index * n_total + i for i in range(n_metal_center)]
            metal_center_positions = extended_lattice.get_positions()[metal_center_indices]

            # Compute distances between ligand atoms and metal center atoms
            distances = np.linalg.norm(
                ligand_positions[:, np.newaxis, :] - metal_center_positions[np.newaxis, :, :],
                axis=2
            )

            # Find the metal center atom with the minimum distance to any ligand atom
            min_distance_idx = np.unravel_index(np.argmin(distances), distances.shape)
            closest_metal_atom_idx = metal_center_indices[min_distance_idx[1]]

            # Remove this metal center atom
            atoms_to_remove.add(closest_metal_atom_idx)

    # Remove the identified atoms and pbc from the extended lattice
    keep_indices = [i for i in range(len(extended_lattice)) if i not in atoms_to_remove]
    extended_lattice = extended_lattice[keep_indices]
    extended_lattice.pbc = False

    ligand_formula = ligand.get_chemical_formula()
    metal_center_formula = metal_center.get_chemical_formula()
    filename = f"{metal_center_formula}_{ligand_formula}_{target_size}_{target_shape}_nanoparticle.cif" # xyz fails for 64 cell and larger
    write(filename, extended_lattice)

    return extended_lattice

