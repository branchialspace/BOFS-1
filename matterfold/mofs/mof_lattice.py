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
    Extends a metal-ligand structure by adding ligands and metal centers to all available atoms
    on the convex hull of every metal center, using the geometric relationships from all
    metal-template pairs in the initial structure.

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

    # Find all metal centers in combined_structure
    metal_indices = []
    current_idx = 0
    while current_idx < len(combined_structure):
        if symbols[current_idx] == metal_symbol:
            metal_indices.append(list(range(current_idx, current_idx + metal_length)))
            current_idx += metal_length
        else:
            current_idx += ligand_length

    if len(metal_indices) < 2:
        raise ValueError("Need at least two metal centers to determine geometric relationships")

    # Get geometric relationships from initial structure
    ligand_start_idx = metal_length
    ligand_positions = positions[ligand_start_idx:ligand_start_idx + ligand_length]

    # Calculate bonding site centroids
    bonding_site1_centroid = np.mean(ligand_positions[bonding_sites[0]], axis=0)
    bonding_site2_centroid = np.mean(ligand_positions[bonding_sites[1]], axis=0)
    ligand_centroid = np.mean(ligand_positions, axis=0)

    # Create extended structure starting with original
    extended_structure = combined_structure.copy()

    # Precompute template structures and their geometric relationships
    templates = []
    template_orientations = []
    for metal_indices_item in metal_indices:
        # Create template by masking current metal center
        template_mask = np.ones(len(combined_structure), dtype=bool)
        template_mask[metal_indices_item] = False
        template_structure = combined_structure[template_mask]
        
        # Calculate metal center centroid
        metal_positions = positions[metal_indices_item]
        metal_centroid = np.mean(metal_positions, axis=0)
        
        # Find template's metal center positions and calculate centroid
        template_metal_indices = []
        template_positions = template_structure.get_positions()
        template_symbols = template_structure.get_chemical_symbols()
        
        current_idx = 0
        while current_idx < len(template_structure):
            if template_symbols[current_idx] == metal_symbol:
                template_metal_indices.append(list(range(current_idx, current_idx + metal_length)))
                current_idx += metal_length
            else:
                current_idx += ligand_length
        
        # Calculate orientations for each template metal center
        orientations = []
        for template_metal_idx in template_metal_indices:
            template_metal_positions = template_positions[template_metal_idx]
            template_metal_centroid = np.mean(template_metal_positions, axis=0)
            
            # Calculate orientation vector from masked metal center to template metal center
            orientation_vector = template_metal_centroid - metal_centroid
            orientations.append({
                'vector': orientation_vector,
                'distance': np.linalg.norm(orientation_vector),
                'template_positions': template_positions
            })
        
        templates.append(template_structure)
        template_orientations.append(orientations)

    # Process each metal center
    for source_metal_idx, source_metal_indices in enumerate(metal_indices):
        # Get source metal center positions
        source_metal_positions = positions[source_metal_indices]
        source_metal_centroid = np.mean(source_metal_positions, axis=0)

        # Get template and its orientations for this metal center
        template_structure = templates[source_metal_idx]
        orientations = template_orientations[source_metal_idx]

        # Get convex hull atoms of source metal center
        hull = ConvexHull(source_metal_positions)
        hull_atoms = np.unique(hull.simplices.flatten())

        # Find the coordinating atoms on the source metal center
        coordinated_atoms = set()
        for other_metal_idx, other_metal_indices in enumerate(metal_indices):
            if other_metal_idx == source_metal_idx:
                continue
            other_centroid = np.mean(positions[other_metal_indices], axis=0)
            distances = np.linalg.norm(source_metal_positions - other_centroid, axis=1)
            coordinated_atoms.add(np.argmin(distances))

        # Remove atoms that are already coordinated
        hull_atoms = np.array([atom for atom in hull_atoms if atom not in coordinated_atoms])

        # For each uncoordinated atom on convex hull
        for hull_atom_idx in hull_atoms:
            hull_atom_pos = source_metal_positions[hull_atom_idx]
            hull_vector = hull_atom_pos - source_metal_centroid

            # For each orientation from the template
            for orientation in orientations:
                # Calculate rotation matrix to align orientation vector with hull vector
                rotation_matrix = calculate_rotation_matrix(orientation['vector'], hull_vector)

                # Create new segment from template
                new_segment = template_structure.copy()
                new_segment_positions = orientation['template_positions'].copy()

                # Transform new segment positions
                closest_coord_atom = min(coordinated_atoms, 
                                      key=lambda x: np.linalg.norm(source_metal_positions[x] - hull_atom_pos))
                translation = hull_atom_pos - source_metal_positions[closest_coord_atom]
                
                new_segment_positions = np.dot(
                    new_segment_positions - source_metal_positions[closest_coord_atom],
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
    """Calculate rotation matrix to align vec1 with vec2 using Rodrigues' rotation formula."""
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
