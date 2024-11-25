# MOF lattice

import numpy as np
from ase import Atoms
from ase.io import write
from scipy.spatial import ConvexHull


def mof_lattice(
    combined_structure: Atoms,
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
) -> Atoms:
    """
    Extends a metal-ligand structure by adding ligands and metal centers to all atoms
    on the convex hull of the first metal center, maintaining the same geometric relationships
    as the initial metal-ligand-metal coordination.

    Parameters:
    - combined_structure: ASE Atoms object of the initial metal-ligand-metal structure
    - ligand: Original ligand Atoms object used to create combined_structure
    - metal_center: Original metal center Atoms object used to create combined_structure
    - bonding_sites: Original list of bonding site indices used in initial docking

    Returns:
    - extended_structure: ASE Atoms object with the extended coordination structure
    """
    # Get positions and symbols
    positions = combined_structure.get_positions()
    symbols = combined_structure.get_chemical_symbols()

    # Identify components in combined_structure
    metal_symbol = metal_center.get_chemical_symbols()[0]
    ligand_length = len(ligand)
    metal_length = len(metal_center)

    # Find metal centers in combined_structure
    metal_indices = []
    current_idx = 0
    while current_idx < len(combined_structure):
        if symbols[current_idx] == metal_symbol:
            metal_indices.append(list(range(current_idx, current_idx + metal_length)))
            current_idx += metal_length
        else:
            current_idx += ligand_length

    extend_from_index = 0  # Always use first metal center
    if len(metal_indices) < 1:
        raise ValueError(f"No metal centers found in structure")

    # Get geometric relationships from initial structure
    ligand_start_idx = metal_length
    ligand_positions = positions[ligand_start_idx:ligand_start_idx + ligand_length]

    # Calculate bonding site centroids
    bonding_site1_centroid = np.mean(ligand_positions[bonding_sites[0]], axis=0)
    bonding_site2_centroid = np.mean(ligand_positions[bonding_sites[1]], axis=0)
    ligand_centroid = np.mean(ligand_positions, axis=0)

    # Get source and target metal centers
    source_metal_indices = metal_indices[extend_from_index]
    source_metal_positions = positions[source_metal_indices]
    source_metal_centroid = np.mean(source_metal_positions, axis=0)

    # Find the coordinating atom on the source metal center (previously find_closest_metal_atom)
    bonding_site_centroid = bonding_site1_centroid if extend_from_index == 0 else bonding_site2_centroid
    distances = np.linalg.norm(source_metal_positions - bonding_site_centroid, axis=1)
    source_coord_idx = np.argmin(distances)

    # Get template structure (everything except source metal center)
    template_mask = np.ones(len(combined_structure), dtype=bool)
    template_mask[source_metal_indices] = False
    template_structure = combined_structure[template_mask]
    template_positions = template_structure.get_positions()

    # Create extended structure starting with original
    extended_structure = combined_structure.copy()

    # Get convex hull atoms of source metal center
    hull = ConvexHull(source_metal_positions)
    hull_atoms = np.unique(hull.simplices.flatten())

    # Remove the atom that's already coordinated
    hull_atoms = hull_atoms[hull_atoms != source_coord_idx]

    # Get vector from source metal centroid to coordinating atom
    source_coord_vector = source_metal_positions[source_coord_idx] - source_metal_centroid

    # For each uncoordinated atom on convex hull
    for hull_atom_idx in hull_atoms:
        hull_atom_pos = source_metal_positions[hull_atom_idx]
        hull_vector = hull_atom_pos - source_metal_centroid

        # Calculate rotation matrix to align hull_vector with source_coord_vector
        rotation_matrix = calculate_rotation_matrix(source_coord_vector, hull_vector)

        # Create new segment from template
        new_segment = template_structure.copy()
        new_segment_positions = template_positions.copy()

        # Transform new segment positions
        translation = hull_atom_pos - source_metal_positions[source_coord_idx]
        new_segment_positions = np.dot(
            new_segment_positions - source_metal_positions[source_coord_idx],
            rotation_matrix.T
        ) + hull_atom_pos

        new_segment.set_positions(new_segment_positions)
        extended_structure += new_segment

    ligand_formula = ligand.get_chemical_formula()
    metal_center_formula = metal_center.get_chemical_formula()
    filename = f"{metal_center_formula}_{ligand_formula}_lattice.xyz"
    write(filename, extended_structure)

    return extended_structure

def calculate_rotation_matrix(vec1, vec2):
    """Calculate rotation matrix to align vec1 with vec2."""
    vec1 = vec1 / np.linalg.norm(vec1)
    vec2 = vec2 / np.linalg.norm(vec2)

    cross_product = np.cross(vec1, vec2)
    dot_product = np.dot(vec1, vec2)

    if np.allclose(dot_product, 1.0):
        return np.eye(3)
    elif np.allclose(dot_product, -1.0):
        # Vectors are antiparallel, rotate 180° around any perpendicular axis
        perpendicular = np.array([1, 0, 0]) if not np.allclose(vec1, [1, 0, 0]) else [0, 1, 0]
        axis = np.cross(vec1, perpendicular)
        axis = axis / np.linalg.norm(axis)
        theta = np.pi
    else:
        axis = cross_product / np.linalg.norm(cross_product)
        theta = np.arccos(dot_product)

    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    rotation_matrix = (np.eye(3) + np.sin(theta) * K +
                      (1 - np.cos(theta)) * np.dot(K, K))

    return rotation_matrix
