# Ligand Metal Docking

import numpy as np
from ase import Atoms
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds


def ligand_metal_docking(
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
    bond_distance: float
) -> Atoms:
    """
    Place metal centers at each bonding site on the ligand using geometric centroids
    and bond angle-aware direction vectors.

    Parameters:
    - ligand: ASE Atoms object of the ligand
    - metal_center: ASE Atoms object of the metal cluster
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site
    - bond_distance: Desired bond distance between the ligand and metal atoms

    Returns:
    - combined_structure: ASE Atoms object with metal centers placed at each bonding site

    Raises:
    - ValueError: If steric hindrance is detected between metal centers or with the ligand
    """
    # Start with a copy of the ligand structure
    combined_structure = ligand.copy()
    ligand_positions = ligand.get_positions()
    previous_metal_positions = []

    # Convert ASE Atoms to RDKit molecule
    symbols = ligand.get_chemical_symbols()
    coords = ligand.get_positions()
    xyz_block = f"{len(symbols)}\n\n"
    for symbol, coord in zip(symbols, coords):
        xyz_block += f"{symbol} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n"

    rdkit_mol = Chem.MolFromXYZBlock(xyz_block)
    Chem.SanitizeMol(rdkit_mol)
    rdDetermineBonds.DetermineConnectivity(rdkit_mol)

    # Loop over each bonding site
    for site_indices in bonding_sites:
        # Get bonding site positions and centroid
        site_positions = ligand_positions[site_indices]
        site_centroid = np.mean(site_positions, axis=0)

        # Calculate direction vector using bond angles
        bond_vectors_sum = np.zeros(3)
        for site_idx in site_indices:
            # Get bonds directly from RDKit molecule
            bonds = rdkit_mol.GetBonds()
            for bond in bonds:
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()

                # If one of the atoms in the bond is our site atom
                if begin_idx == site_idx or end_idx == site_idx:
                    neighbor_idx = end_idx if begin_idx == site_idx else begin_idx
                    if neighbor_idx not in site_indices:
                        # Calculate bond vector from bonding site atom to neighbor
                        bond_vector = ligand_positions[site_idx] - ligand_positions[neighbor_idx]
                        bond_vector = bond_vector / np.linalg.norm(bond_vector)
                        bond_vectors_sum += bond_vector

        # The direction vector is the sum of normalized bond vectors
        direction_vector = bond_vectors_sum
        direction_vector = direction_vector / np.linalg.norm(direction_vector)

        # Define objective function to find optimal bonding position
        def position_objective_function(scale_factor):
            proposed_position = site_centroid + scale_factor * direction_vector
            # Calculate distances from proposed_position to each atom in the bonding site
            distances = np.linalg.norm(site_positions - proposed_position, axis=1)
            # Objective: minimize the sum of squared deviations from bond_distance
            return np.sum((distances - bond_distance) ** 2)

        # Optimize the scaling factor to find optimal bonding position
        position_result = minimize(position_objective_function, bond_distance, method='L-BFGS-B')
        optimal_position = site_centroid + position_result.x[0] * direction_vector

        # Identify coordinating atom from metal center
        metal_positions = metal_center.get_positions()
        hull = ConvexHull(metal_positions)
        coordinating_atom_index = np.unique(hull.simplices.flatten())[0]

        # Create a copy of metal center and align it
        metal_center_copy = metal_center.copy()
        metal_center_positions = metal_center_copy.get_positions()

        # Translate metal center to align coordinating atom with optimal position
        translation_vector = optimal_position - metal_center_positions[coordinating_atom_index]
        metal_center_positions += translation_vector

        # Optimize rotation around the binding position
        def rotation_objective_function(rotation_angles):
            theta_x, theta_y, theta_z = rotation_angles
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

            rotated_positions = np.dot(
                metal_center_positions - optimal_position,
                R.T
            ) + optimal_position

            # Calculate distances to ligand and previous metal centers
            metal_indices = np.arange(len(metal_center))
            metal_indices = np.delete(metal_indices, coordinating_atom_index)
            metal_atoms_positions = rotated_positions[metal_indices]

            # Distances to ligand
            distances_ligand = np.linalg.norm(
                metal_atoms_positions[:, np.newaxis, :] -
                ligand_positions[np.newaxis, :, :],
                axis=2
            )
            min_distances = np.min(distances_ligand, axis=1)

            # Distances to previous metal centers
            if previous_metal_positions:
                prev_metals = np.vstack(previous_metal_positions)
                distances_metals = np.linalg.norm(
                    metal_atoms_positions[:, np.newaxis, :] -
                    prev_metals[np.newaxis, :, :],
                    axis=2
                )
                min_distances = np.minimum(
                    min_distances,
                    np.min(distances_metals, axis=1)
                )

            # We want to maximize the minimum distance, so we return the negative sum
            return -np.sum(min_distances)

        # Optimize rotation
        initial_angles = np.array([0.0, 0.0, 0.0])
        rotation_result = minimize(
            rotation_objective_function,
            initial_angles,
            method='L-BFGS-B'
        )

        # Apply optimal rotation
        theta_x, theta_y, theta_z = rotation_result.x
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

        final_positions = np.dot(
            metal_center_positions - optimal_position,
            R.T
        ) + optimal_position

        # Check for steric hindrance
        metal_indices = np.arange(len(metal_center))
        metal_indices = np.delete(metal_indices, coordinating_atom_index)
        metal_atoms_positions = final_positions[metal_indices]

        distances_ligand = np.linalg.norm(
            metal_atoms_positions[:, np.newaxis, :] -
            ligand_positions[np.newaxis, :, :],
            axis=2
        )
        min_distances = np.min(distances_ligand, axis=1)

        if previous_metal_positions:
            prev_metals = np.vstack(previous_metal_positions)
            distances_metals = np.linalg.norm(
                metal_atoms_positions[:, np.newaxis, :] -
                prev_metals[np.newaxis, :, :],
                axis=2
            )
            min_distances = np.minimum(
                min_distances,
                np.min(distances_metals, axis=1)
            )

        if np.any(min_distances < bond_distance):
            raise ValueError("Steric hindrance detected between metal centers or with the ligand.")

        # Set final positions and update structure
        metal_center_copy.set_positions(final_positions)
        combined_structure += metal_center_copy
        previous_metal_positions.append(metal_atoms_positions)

    return combined_structure
        previous_metal_positions.append(metal_atoms_positions)

    return combined_structure
