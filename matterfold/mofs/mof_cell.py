# MOF primitive cubic/ rhombohedral lattice unit cell

import numpy as np
from ase import Atoms
from ase.io import write
from scipy.spatial import ConvexHull


def mof_cell(
    combined_structure: Atoms,
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list
) -> Atoms:
    """
    Constructs a unit cell by replicating the ligand and metal centers to fully coordinate all the
    convex hull atoms of the metal centers with ligands, forming a primitive cubic/ rhombohedral lattice.
    Periodic boundary conditions are applied to maintain proper bonding across unit cells.

    Parameters:
    - combined_structure: ASE Atoms object returned by ligand_metal_docking (ligand + two metal centers).
    - ligand: ASE Atoms object of the ligand.
    - metal_center: ASE Atoms object of the metal cluster.
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site on the ligand.

    Returns:
    - unit_cell_structure: ASE Atoms object representing the unit cell with PBC.
    """
    # Get the number of atoms in the ligand and metal_center
    num_ligand_atoms = len(ligand)
    num_metal_atoms = len(metal_center)

    # Total expected number of atoms: ligand + 2 * metal_center
    expected_total_atoms = num_ligand_atoms + 2 * num_metal_atoms

    if len(combined_structure) != expected_total_atoms:
        raise ValueError(f"Expected total atoms: {expected_total_atoms}, but got {len(combined_structure)}")

    # Identify indices of ligand and metal atoms in combined_structure
    ligand_indices = range(num_ligand_atoms)
    metal1_indices = range(num_ligand_atoms, num_ligand_atoms + num_metal_atoms)
    metal2_indices = range(num_ligand_atoms + num_metal_atoms, expected_total_atoms)

    # Get positions and symbols
    positions = combined_structure.get_positions()
    symbols = combined_structure.get_chemical_symbols()

    # Get positions of the two metal centers
    metal1_positions = positions[metal1_indices]
    metal2_positions = positions[metal2_indices]

    # Calculate the centroids of the metal centers
    metal1_centroid = np.mean(metal1_positions, axis=0)
    metal2_centroid = np.mean(metal2_positions, axis=0)

    # Calculate the vector between the two metal centers
    vector_metal_centers = metal2_centroid - metal1_centroid
    distance_metal_centers = np.linalg.norm(vector_metal_centers)
    lattice_constant = distance_metal_centers # Define lattice constant

    # Define lattice vectors for a cubic lattice
    a1 = lattice_constant * np.array([1, 0, 0])
    a2 = lattice_constant * np.array([0, 1, 0])
    a3 = lattice_constant * np.array([0, 0, 1])
    unit_cell = np.array([a1, a2, a3])

    # Create an empty Atoms object for the unit cell
    unit_cell_structure = Atoms(cell=unit_cell, pbc=[True, True, True])

    # Place the metal center at the origin
    metal_center_positions = metal_center.get_positions()
    metal_center_symbols = metal_center.get_chemical_symbols()
    metal_center_centroid = np.mean(metal_center_positions, axis=0)
    metal_center_positions -= metal_center_centroid  # Center at origin

    # Identify atoms on the convex hull of the metal center
    hull = ConvexHull(metal_center_positions)
    convex_hull_indices = np.unique(hull.simplices.flatten())

    # Place ligands around the metal center at each bonding site
    for idx in convex_hull_indices:
        site_position = metal_center_positions[idx]
        site_vector = site_position - metal_center_centroid

        # Rotate the ligand to align its bonding site with the site_vector
        ligand_copy = ligand.copy()

        # Compute the ligand centroid
        ligand_centroid = np.mean(ligand_copy.get_positions(), axis=0)
        # Compute the bonding site centroid
        bonding_site_positions = ligand_copy.get_positions()[bonding_sites[0]]
        bonding_site_centroid = np.mean(bonding_site_positions, axis=0)
        # Compute the vector from ligand centroid to bonding site centroid
        ligand_axis = bonding_site_centroid - ligand_centroid
        ligand_axis /= np.linalg.norm(ligand_axis)

        # Normalize target_vector
        target_vector = site_vector.copy()
        target_vector /= np.linalg.norm(target_vector)

        # Calculate the rotation axis and angle
        rotation_axis = np.cross(ligand_axis, target_vector)
        sin_angle = np.linalg.norm(rotation_axis)
        cos_angle = np.dot(ligand_axis, target_vector)
        rotation_angle = np.arctan2(sin_angle, cos_angle)

        if sin_angle >= 1e-6:
            rotation_axis /= sin_angle  # Normalize rotation axis
            cos_angle = np.cos(rotation_angle)
            sin_angle = np.sin(rotation_angle)
            ux, uy, uz = rotation_axis

            rotation_matrix = np.array([
                [cos_angle + ux**2 * (1 - cos_angle),
                 ux * uy * (1 - cos_angle) - uz * sin_angle,
                 ux * uz * (1 - cos_angle) + uy * sin_angle],
                [uy * ux * (1 - cos_angle) + uz * sin_angle,
                 cos_angle + uy**2 * (1 - cos_angle),
                 uy * uz * (1 - cos_angle) - ux * sin_angle],
                [uz * ux * (1 - cos_angle) - uy * sin_angle,
                 uz * uy * (1 - cos_angle) + ux * sin_angle,
                 cos_angle + uz**2 * (1 - cos_angle)]
            ])

            # Apply rotation around the ligand centroid
            positions = ligand_copy.get_positions() - ligand_centroid
            rotated_positions = positions.dot(rotation_matrix.T)
            ligand_copy.set_positions(rotated_positions + ligand_centroid)

        # Translate the ligand so that its bonding site aligns with the convex hull atom position
        bonding_site_positions = ligand_copy.get_positions()[bonding_sites[0]]
        bonding_site_centroid = np.mean(bonding_site_positions, axis=0)
        translation = site_position - bonding_site_centroid
        ligand_copy.translate(translation)

        # Add the rotated and translated ligand to the unit cell structure
        unit_cell_structure += ligand_copy

    # Add the metal center to the unit cell structure
    unit_cell_structure.extend(Atoms(
        symbols=metal_center_symbols,
        positions=metal_center_positions,
    ))

    # Apply periodic boundary conditions
    unit_cell_structure.wrap()

    ligand_formula = ligand.get_chemical_formula()
    metal_center_formula = metal_center.get_chemical_formula()
    filename = f"{metal_center_formula}_{ligand_formula}_cell.xyz"
    write(filename, unit_cell_structure)


    return unit_cell_structure
