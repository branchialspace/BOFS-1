# Ligand Metal Docking

import numpy as np
from ase import Atoms
from scipy.spatial import ConvexHull
from scipy.optimize import minimize


def ligand_metal_docking(
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
    bond_distance: float
) -> Atoms:
    """
    Place metal centers at each bonding site on the ligand using the Kabsch algorithm and a force-based rotation optimization.

    Parameters:
    - ligand: ASE Atoms object of the ligand.
    - metal_center: ASE Atoms object of the metal cluster.
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site.
    - bond_distance: Desired bond distance between the ligand and metal atoms.
    - coordinating_atom_index: Optional; index of the atom in the metal center to act as the coordinating atom.
                               If None, an atom from the convex hull will be chosen arbitrarily.

    Returns:
    - combined_structure: ASE Atoms object with metal centers placed at each bonding site.

    Raises:
    - ValueError: If steric hindrance is detected between metal centers or with the ligand.
    """
    # Start with a copy of the ligand structure
    combined_structure = ligand.copy()
    ligand_positions = ligand.get_positions()
    ligand_indices = np.arange(len(ligand))

    # Initialize a list to keep track of metal center positions (excluding coordinating atoms)
    previous_metal_positions = []

    # Loop over each bonding site
    for site_indices in bonding_sites:
        # Extract positions of the ligand's bonding-site atoms
        bonding_positions = ligand_positions[site_indices]

        # Compute centroid of the bonding-site atoms
        ligand_centroid = np.mean(bonding_positions, axis=0)

        # Compute normal vector to the plane defined by the bonding-site atoms (if at least 3 atoms)
        if len(bonding_positions) >= 3:
            vec1 = bonding_positions[1] - bonding_positions[0]
            vec2 = bonding_positions[2] - bonding_positions[0]
            normal_vector = np.cross(vec1, vec2)
            normal_vector /= np.linalg.norm(normal_vector)
        else:
            # If less than 3 atoms, define an arbitrary normal vector
            normal_vector = np.array([0.0, 0.0, 1.0])

        # Identify the coordinating atom
        metal_positions = metal_center.get_positions()
        # Compute the convex hull of the metal center to find surface atoms
        hull = ConvexHull(metal_positions)
        convex_hull_indices = np.unique(hull.simplices.flatten())
        # Choose any atom from the convex hull as the coordinating atom
        coordinating_atom_index = convex_hull_indices[0]

        coordinating_atom_position = metal_positions[coordinating_atom_index]

        # Create a copy of the metal center and translate it so that the coordinating atom is at the origin
        metal_center_copy = metal_center.copy()
        metal_center_positions = metal_center_copy.get_positions()
        metal_center_positions -= coordinating_atom_position
        metal_center_copy.set_positions(metal_center_positions)

        # Define the desired position of the coordinating atom
        # Optimize the position C to minimize the sum of squared differences between distances and bond_distance
        def objective_function(C):
            distances = np.linalg.norm(bonding_positions - C, axis=1)
            return np.sum((distances - bond_distance) ** 2)

        # Initial guess for C is at the ligand centroid shifted along the normal vector
        initial_guess = ligand_centroid + bond_distance * normal_vector

        # Optimize the position of the coordinating atom
        result = minimize(objective_function, initial_guess, method='L-BFGS-B')
        coordinating_atom_new_position = result.x

        # Place the coordinating atom at the optimized position
        translation_vector = coordinating_atom_new_position

        # Apply the translation to the metal center positions
        metal_center_positions += translation_vector
        metal_center_copy.set_positions(metal_center_positions)

        # Rotate the metal center around the coordinating atom to maximize the distance from the ligand and previous metal centers
        # Define a force-based function to optimize the rotation
        def rotation_objective_function(rotation_angles):
            # rotation_angles = [theta_x, theta_y, theta_z]
            theta_x, theta_y, theta_z = rotation_angles
            # Create rotation matrices around x, y, z axes
            Rx = np.array([[1, 0, 0],
                           [0, np.cos(theta_x), -np.sin(theta_x)],
                           [0, np.sin(theta_x), np.cos(theta_x)]])
            Ry = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
                           [0, 1, 0],
                           [-np.sin(theta_y), 0, np.cos(theta_y)]])
            Rz = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
                           [np.sin(theta_z), np.cos(theta_z), 0],
                           [0, 0, 1]])
            # Combined rotation matrix
            R = Rz @ Ry @ Rx
            # Rotate the metal center positions around the coordinating atom (which is at coordinating_atom_new_position)
            rotated_positions = np.dot(metal_center_positions - coordinating_atom_new_position, R.T) + coordinating_atom_new_position
            # Exclude the coordinating atom from steric calculations
            metal_indices = np.arange(len(metal_center))
            metal_indices = np.delete(metal_indices, coordinating_atom_index)
            metal_atoms_positions = rotated_positions[metal_indices]
            # Compute distances between metal atoms and ligand atoms
            distances_ligand = np.linalg.norm(metal_atoms_positions[:, np.newaxis, :] - ligand_positions[np.newaxis, :, :], axis=2)
            min_distances_ligand = np.min(distances_ligand, axis=1)
            # Initialize total minimal distances with distances to ligand
            min_distances = min_distances_ligand.copy()
            # If there are previous metal centers, compute distances to them
            if previous_metal_positions:
                previous_metal_positions_array = np.vstack(previous_metal_positions)
                distances_previous_metals = np.linalg.norm(metal_atoms_positions[:, np.newaxis, :] - previous_metal_positions_array[np.newaxis, :, :], axis=2)
                min_distances_previous_metals = np.min(distances_previous_metals, axis=1)
                # Update min_distances to include distances to previous metals
                min_distances = np.minimum(min_distances, min_distances_previous_metals)
            # Objective is to maximize the minimal distances (minimize the negative sum)
            return -np.sum(min_distances)
        
        # Initial rotation angles (no rotation)
        initial_rotation_angles = np.array([0.0, 0.0, 0.0])

        # Optimize the rotation angles
        rotation_result = minimize(rotation_objective_function, initial_rotation_angles, method='L-BFGS-B')
        optimized_angles = rotation_result.x
        theta_x, theta_y, theta_z = optimized_angles

        # Compute the optimized rotation matrix
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(theta_x), -np.sin(theta_x)],
                       [0, np.sin(theta_x), np.cos(theta_x)]])
        Ry = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
                       [0, 1, 0],
                       [-np.sin(theta_y), 0, np.cos(theta_y)]])
        Rz = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
                       [np.sin(theta_z), np.cos(theta_z), 0],
                       [0, 0, 1]])
        R = Rz @ Ry @ Rx

        # Apply the optimized rotation to the metal center positions
        rotated_positions = np.dot(metal_center_positions - coordinating_atom_new_position, R.T) + coordinating_atom_new_position
        metal_center_copy.set_positions(rotated_positions)

        # Now check for steric hindrance
        # Exclude the coordinating atom from steric calculations
        metal_indices = np.arange(len(metal_center))
        metal_indices = np.delete(metal_indices, coordinating_atom_index)
        metal_atoms_positions = rotated_positions[metal_indices]

        # Compute distances between metal atoms and ligand atoms
        distances_ligand = np.linalg.norm(metal_atoms_positions[:, np.newaxis, :] - ligand_positions[np.newaxis, :, :], axis=2)
        min_distances_ligand = np.min(distances_ligand, axis=1)

        # Initialize total minimal distances with distances to ligand
        min_distances = min_distances_ligand.copy()

        # If there are previous metal centers, compute distances to them
        if previous_metal_positions:
            previous_metal_positions_array = np.vstack(previous_metal_positions)
            distances_previous_metals = np.linalg.norm(metal_atoms_positions[:, np.newaxis, :] - previous_metal_positions_array[np.newaxis, :, :], axis=2)
            min_distances_previous_metals = np.min(distances_previous_metals, axis=1)
            # Update min_distances to include distances to previous metals
            min_distances = np.minimum(min_distances, min_distances_previous_metals)

        # Define a threshold for steric hindrance, e.g., 1.0 Å (adjust as needed)
        steric_threshold = 1.0  # Angstroms

        # Check if any minimal distance is below the threshold
        if np.any(min_distances < steric_threshold):
            raise ValueError("Steric hindrance detected between metal centers or with the ligand.")

        # Combine the metal center with the combined structure
        combined_structure += metal_center_copy

        # Add the positions of the current metal center (excluding coordinating atom) to previous_metal_positions
        previous_metal_positions.append(metal_atoms_positions)

    return combined_structure
