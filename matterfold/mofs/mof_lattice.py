# MOF lattice

import numpy as np
from ase import Atoms
from ase.io import write
from scipy.spatial import ConvexHull


def extend_metal_structure(
    combined_structure: Atoms,
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
) -> Atoms:
    """
    Extends a metal-ligand structure by adding ligands and metal centers to all atoms
    on the convex hull of a selected metal center, maintaining the same geometric relationships
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
            
    # Identify which metal center atoms are coordinated to which bonding sites
    def find_closest_metal_atom(bonding_site_centroid, metal_pos):
        distances = np.linalg.norm(metal_pos - bonding_site_centroid, axis=1)
        return np.argmin(distances)
    
    # Get geometric relationships from initial structure
    ligand_start_idx = metal_length
    ligand_positions = positions[ligand_start_idx:ligand_start_idx + ligand_length]
    
    # Calculate bonding site centroids
    bonding_site1_centroid = np.mean(ligand_positions[bonding_sites[0]], axis=0)
    bonding_site2_centroid = np.mean(ligand_positions[bonding_sites[1]], axis=0)
    ligand_centroid = np.mean(ligand_positions, axis=0)
    
    # Find metal centers' positions and coordinating atoms
    metal1_positions = positions[metal_indices[0]]
    metal1_centroid = np.mean(metal1_positions, axis=0)
    metal1_coord_idx = find_closest_metal_atom(bonding_site1_centroid, metal1_positions)
    
    metal2_positions = positions[metal_indices[1]]
    metal2_centroid = np.mean(metal2_positions, axis=0)
    metal2_coord_idx = find_closest_metal_atom(bonding_site2_centroid, metal2_positions)
    
    # Calculate relative vectors and angles for metal1
    metal1_coord_vector = metal1_positions[metal1_coord_idx] - metal1_centroid
    site1_vector = bonding_site1_centroid - metal1_positions[metal1_coord_idx]
    
    # Create extended structure starting with original
    extended_structure = combined_structure.copy()
    
    # For each metal center in the structure
    for metal_center_indices in metal_indices:
        metal_positions = positions[metal_center_indices]
        metal_centroid = np.mean(metal_positions, axis=0)
        
        # Get convex hull atoms
        hull = ConvexHull(metal_positions)
        hull_atoms = np.unique(hull.simplices.flatten())
        
        # Remove the atom that's already coordinated
        coord_atom = find_closest_metal_atom(bonding_site1_centroid, metal_positions)
        hull_atoms = hull_atoms[hull_atoms != coord_atom]
        
        # For each atom on convex hull
        for hull_atom_idx in hull_atoms:
            # Calculate new structure position based on relative geometry
            hull_atom_pos = metal_positions[hull_atom_idx]
            hull_vector = hull_atom_pos - metal_centroid
            
            # Calculate rotation matrix to align hull_vector with original metal1_coord_vector
            rotation_matrix = calculate_rotation_matrix(metal1_coord_vector, hull_vector)
            
            # Create new structure segment
            new_segment = combined_structure.copy()
            new_segment_positions = new_segment.get_positions()
            
            # Remove the metal center we're extending from
            mask = np.ones(len(new_segment), dtype=bool)
            mask[metal_center_indices] = False
            new_segment = new_segment[mask]
            new_segment_positions = new_segment_positions[mask]
            
            # Translate and rotate new segment
            translation = hull_atom_pos - metal1_positions[metal1_coord_idx]
            new_segment_positions = np.dot(
                new_segment_positions - metal1_positions[metal1_coord_idx],
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
