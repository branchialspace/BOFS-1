# Bisected bonding-site charge symmetry for 1d bridging ligands
# Steric hindrance in single-cluster octahedral coordination
    # Bonding site angle = max distance from internal ligand bonds

import numpy as np
from ase import Atoms
from ase.io import write

def find_bonding_site_symmetry_axis(ligand: Atoms, bonding_sites: dict):
    """
    Identify the axis of symmetry of bonding sites in bridging ligands by bisecting the ligand
    to equalize the density of unbonded electrons on either side while maximizing the distance
    between charge densities.

    Parameters:
    - ligand: ASE Atoms object representing the ligand molecule.
    - bonding_sites: Dictionary containing information about each atom's bonding role.

    Returns:
    - best_plane: Dictionary containing the normal vector and point defining the bisecting plane.
    """
    # Extract positions and unbonded electron counts
    positions = ligand.get_positions()
    num_atoms = len(positions)
    weights = np.array([bonding_sites[i]['non_bonded_electrons'] for i in range(num_atoms)])
    weights = np.where(weights > 0, weights, 0)

    # Normalize weights to avoid division by zero
    total_weight = np.sum(weights)
    if total_weight == 0:
        raise ValueError("Total unbonded electrons are zero after disregarding negatives.")

    # Compute the weighted center of unbonded electrons
    weighted_positions = positions * weights[:, np.newaxis]
    center = np.sum(weighted_positions, axis=0) / total_weight

    # Center the positions
    centered_positions = positions - center

    # Perform PCA on the weighted positions
    # Compute the covariance matrix manually since np.cov can't handle negative weights
    deviations = centered_positions
    covariance_matrix = np.zeros((3, 3))
    for i in range(num_atoms):
        if weights[i] > 0:
            deviation = deviations[i][:, np.newaxis]
            covariance_matrix += weights[i] * deviation @ deviation.T
    covariance_matrix /= total_weight

    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

    # The principal component is the eigenvector with the largest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # The normal vector of the bisecting plane is orthogonal to the principal axis
    # We can take the eigenvector corresponding to the smallest eigenvalue as the normal vector
    normal_vector = eigenvectors[:, np.argmin(eigenvalues)]
    normal_vector /= np.linalg.norm(normal_vector)

    # Project positions onto the normal vector
    projections = centered_positions @ normal_vector

    # Sort projections and corresponding weights
    sorted_indices = np.argsort(projections)
    sorted_projections = projections[sorted_indices]
    sorted_weights = weights[sorted_indices]
    cumulative_weights = np.cumsum(sorted_weights)

    # Find the plane position that balances the weights
    idx = np.searchsorted(cumulative_weights, total_weight / 2)
    if idx == 0 or idx == len(projections):
        d = projections[idx]
    else:
        # Average the two projections
        d = (projections[idx - 1] + projections[idx]) / 2

    # The point on the plane is center + normal_vector * d
    point_on_plane = center + normal_vector * d

    # Return the plane parameters
    best_plane = {
        'normal_vector': normal_vector,
        'point_on_plane': point_on_plane
    }

    return best_plane

def visualize_symmetry_axis(ligand: Atoms, symmetry_axis: dict):
    """
    Visualize the symmetry axis and ligand by creating a modified XYZ file for Avogadro2.

    Parameters:
    - ligand: ASE Atoms object representing the ligand molecule.
    - symmetry_axis: Dictionary containing 'normal_vector' and 'point_on_plane'.
    """
    # Extract normal vector and a point on the plane
    n = symmetry_axis['normal_vector']
    point = symmetry_axis['point_on_plane']

    # Calculate two points along the normal vector to represent the axis
    # Extend the axis beyond the ligand for better visualization
    length = np.linalg.norm(ligand.get_positions().ptp(axis=0)) * 1.5  # 1.5 times the ligand size
    axis_point1 = point + n * length
    axis_point2 = point - n * length

    # Create dummy atoms at these points
    # Using 'X' as the element symbol for dummy atoms (ASE recognizes 'X' with atomic number 0)
    dummy_symbols = ['X', 'X']
    dummy_positions = [axis_point1, axis_point2]
    dummy_atoms = Atoms(symbols=dummy_symbols, positions=dummy_positions)

    # Combine ligand atoms with dummy atoms
    combined_atoms = ligand + dummy_atoms

    # Set the atomic numbers to zero to distinguish them
    for atom in combined_atoms[-2:]:
        atom.number = 0

    # Write the combined atoms to an XYZ file
    filename = "".join(ligand.get_chemical_symbols()) + "_bisected.xyz"
    write(filename, combined_atoms)

    print(f"Visualization file '{filename}' has been created.")
