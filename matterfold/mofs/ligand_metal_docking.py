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
        
        # Step 2: Identify the metal atoms that will coordinate
        # Compute the convex hull of the metal cluster to find surface atoms
        metal_positions = metal_center_copy.get_positions()
        hull = ConvexHull(metal_positions)
        surface_atom_indices = np.unique(hull.simplices.flatten())
        surface_atom_positions = metal_positions[surface_atom_indices]
        
        # Compute the centroid of the ligand
        ligand_centroid = ligand.get_center_of_mass()
        
        # Vector from ligand centroid to bonding-site atoms
        bonding_vectors = bonding_positions - ligand_centroid
        if len(bonding_vectors) > 1:
            # Use the average direction if multiple atoms
            bonding_direction = np.mean(bonding_vectors, axis=0)
        else:
            bonding_direction = bonding_vectors[0]
        
        # Normalize the bonding direction
        norm = np.linalg.norm(bonding_direction)
        if norm == 0:
            raise ValueError("Bonding direction vector has zero magnitude.")
        bonding_direction /= norm
        
        # Vector from metal centroid to each surface atom
        metal_centroid = np.mean(metal_positions, axis=0)
        metal_vectors = surface_atom_positions - metal_centroid
        metal_vectors_normalized = metal_vectors / np.linalg.norm(metal_vectors, axis=1)[:, np.newaxis]
        
        # Compute the angle between bonding direction and each metal vector
        angles = np.arccos(np.clip(np.dot(metal_vectors_normalized, bonding_direction), -1.0, 1.0))
        
        # Select the metal atom whose vector is most aligned with the bonding direction
        coordinating_atom_index = surface_atom_indices[np.argmin(angles)]
        coordinating_position = metal_positions[coordinating_atom_index]
        
        # Step 3: Compute the optimal rotation and translation using the Kabsch algorithm
        # Number of bonding-site atoms
        N = len(bonding_positions)
        
        # Select N surface atoms on the metal cluster
        if N == 1:
            # If only one bonding-site atom, use the coordinating atom
            metal_indices = [coordinating_atom_index]
            metal_positions_subset = metal_positions[metal_indices]
        else:
            # Compute distances from coordinating atom to surface atoms
            distances = np.linalg.norm(surface_atom_positions - coordinating_position, axis=1)
            # Include the coordinating atom itself
            closest_indices = np.argsort(distances)[:N]
            metal_indices = surface_atom_indices[closest_indices]
            metal_positions_subset = metal_positions[metal_indices]
        
        # Compute centroids
        ligand_centroid_subset = np.mean(bonding_positions, axis=0)
        metal_centroid_subset = np.mean(metal_positions_subset, axis=0)
        
        # Center the positions
        P = bonding_positions - ligand_centroid_subset
        Q = metal_positions_subset - metal_centroid_subset
        
        # Compute covariance matrix
        C = np.dot(Q.T, P)
        
        # Compute SVD
        U, S, Vt = np.linalg.svd(C)
        
        # Compute rotation matrix
        D = np.diag([1, 1, np.linalg.det(np.dot(Vt.T, U.T))])
        R = np.dot(np.dot(Vt.T, D), U.T)
        
        # Apply rotation to the entire metal center
        rotated_metal_positions = np.dot(metal_positions - metal_centroid, R)
        
        # Calculate the translation vector
        # Adjust the distance between the centroids along the bonding direction
        translation = ligand_centroid_subset + bond_distance * bonding_direction - (np.dot(metal_centroid, R))
        
        # Apply translation
        new_metal_positions = rotated_metal_positions + translation
        
        # Update the metal center positions
        metal_center_copy.set_positions(new_metal_positions)
        
        # Step 4: Combine the metal center with the combined structure
        combined_structure += metal_center_copy
        
    return combined_structure
