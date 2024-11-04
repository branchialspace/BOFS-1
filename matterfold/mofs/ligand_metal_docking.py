# Ligand Metal Docking

import numpy as np
from ase import Atoms


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
        # For simplicity, select the closest metal atoms to the bonding site(s)
        metal_positions = metal_center_copy.get_positions()
        distances = np.linalg.norm(
            metal_positions[:, np.newaxis, :] - bonding_positions[np.newaxis, :, :], axis=2
        )
        coordinating_atom_indices = np.argmin(distances, axis=0)
        coordinating_positions = metal_positions[coordinating_atom_indices]
        
        # Step 3: Compute the optimal rotation and translation using the Kabsch algorithm
        # Center the coordinates
        ligand_centroid = np.mean(bonding_positions, axis=0)
        metal_centroid = np.mean(coordinating_positions, axis=0)
        ligand_coords_centered = bonding_positions - ligand_centroid
        metal_coords_centered = coordinating_positions - metal_centroid
        
        # Compute covariance matrix
        covariance_matrix = np.dot(metal_coords_centered.T, ligand_coords_centered)
        
        # Singular Value Decomposition
        V, S, Wt = np.linalg.svd(covariance_matrix)
        d = np.linalg.det(np.dot(V, Wt))
        D = np.diag([1, 1, d])
        
        # Rotation matrix
        rotation_matrix = np.dot(np.dot(V, D), Wt)
        
        # Apply rotation to the entire metal center
        rotated_metal_positions = np.dot(metal_center_copy.get_positions() - metal_centroid, rotation_matrix)
        
        # Step 4: Translate the metal center to the bonding site
        # Adjust the translation to match the desired bond distance
        current_distance = np.linalg.norm(
            (rotated_metal_positions[coordinating_atom_indices] + metal_centroid) - bonding_positions, axis=1
        ).mean()
        distance_correction = bond_distance - current_distance
        # Direction vector from metal centroid to ligand centroid
        direction_vector = ligand_centroid - metal_centroid
        # Normalize direction vector
        if np.linalg.norm(direction_vector) > 0:
            direction_vector /= np.linalg.norm(direction_vector)
        else:
            direction_vector = np.zeros(3)
        # Apply distance correction along the direction vector
        translation_vector = ligand_centroid - metal_centroid + distance_correction * direction_vector
        
        # Apply rotation and translation
        new_metal_positions = rotated_metal_positions + metal_centroid + translation_vector
        
        # Update the metal center positions
        metal_center_copy.set_positions(new_metal_positions)
        
        # Step 5: Combine the metal center with the combined structure
        combined_structure += metal_center_copy
    
    return combined_structure
