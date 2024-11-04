# Ligand Metal Docking

import numpy as np
from ase import Atoms
from scipy.spatial import ConvexHull


def ligand_metal_docking(
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
    bond_distance: float
) -> Atoms:
    """
    Place metal centers at each bonding site cluster on the ligand using the Kabsch algorithm RMSD method.
    
    Parameters:
    - ligand: ASE Atoms object of the ligand.
    - metal_center: ASE Atoms object of the metal cluster.
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site.
    - bond_distance: Desired bond distance between the ligand and metal atoms.
    
    Returns:
    - combined_structure: ASE Atoms object with metal centers placed at each bonding site.
    """
    # Start with the ligand structure
    combined_structure = ligand.copy()
    
    # Loop over each bonding site cluster
    for bonding_site_atoms in bonding_sites:
        # Copy the metal center for placement
        metal_center_copy = metal_center.copy()
        
        # Step 1: Extract positions of the ligand's bonding-site atoms
        ligand_positions = ligand.get_positions()
        bonding_positions = ligand_positions[bonding_site_atoms]
        
        # Step 2: Identify the metal atom that will coordinate
        # Compute the convex hull of the metal cluster to find surface atoms
        metal_positions = metal_center_copy.get_positions()
        hull = ConvexHull(metal_positions)
        surface_atom_indices = np.unique(hull.simplices.flatten())
        surface_atom_positions = metal_positions[surface_atom_indices]
        
        # For simplicity, select the surface atom that best aligns along the computed angle
        # Compute the centroid of the bonding-site atoms
        ligand_centroid = np.mean(bonding_positions, axis=0)
        
        # Vector from ligand centroid to bonding-site centroid
        bonding_vector = bonding_positions - ligand_centroid
        if len(bonding_vector) > 1:
            # Use the average direction if multiple atoms
            bonding_direction = np.mean(bonding_vector, axis=0)
        else:
            bonding_direction = bonding_vector[0]
        bonding_direction /= np.linalg.norm(bonding_direction)
        
        # Vector from metal centroid to each surface atom
        metal_centroid = np.mean(metal_positions, axis=0)
        metal_vectors = surface_atom_positions - metal_centroid
        metal_vectors_normalized = metal_vectors / np.linalg.norm(metal_vectors, axis=1)[:, np.newaxis]
        
        # Compute the angle between bonding direction and each metal vector
        angles = np.arccos(np.clip(np.dot(metal_vectors_normalized, bonding_direction), -1.0, 1.0))
        
        # Select the metal atom whose vector is most aligned with the bonding direction
        coordinating_atom_index = surface_atom_indices[np.argmin(angles)]
        coordinating_position = metal_positions[coordinating_atom_index]
        
        # Step 3: Compute the optimal rotation and translation
        # We need at least three non-colinear points to compute a rotation
        # Since we have only one coordinating atom, we'll align the vectors
        # Vector from metal centroid to coordinating atom
        metal_vector = coordinating_position - metal_centroid
        metal_vector /= np.linalg.norm(metal_vector)
        
        # Compute rotation matrix to align metal_vector to bonding_direction
        v = np.cross(metal_vector, bonding_direction)
        c = np.dot(metal_vector, bonding_direction)
        s = np.linalg.norm(v)
        if s != 0:
            kmat = np.array([[0, -v[2], v[1]],
                             [v[2], 0, -v[0]],
                             [-v[1], v[0], 0]])
            rotation_matrix = np.eye(3) + kmat + kmat @ kmat * ((1 - c) / (s ** 2))
        else:
            rotation_matrix = np.eye(3)  # Vectors are aligned
        
        # Apply rotation to the entire metal center
        rotated_metal_positions = np.dot(metal_positions - metal_centroid, rotation_matrix)
        
        # Step 4: Translate the metal center to the bonding site
        # Place the coordinating atom at the desired bond distance along the bonding direction
        translation_vector = ligand_centroid + bond_distance * bonding_direction - (metal_centroid + np.dot(coordinating_position - metal_centroid, rotation_matrix))
        
        # Apply rotation and translation
        new_metal_positions = rotated_metal_positions + metal_centroid + translation_vector
        
        # Update the metal center positions
        metal_center_copy.set_positions(new_metal_positions)
        
        # Step 5: Combine the metal center with the combined structure
        combined_structure += metal_center_copy
        
    return combined_structure
