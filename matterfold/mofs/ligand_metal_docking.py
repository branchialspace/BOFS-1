# Ligand Metal Docking

import numpy as np
from ase import Atoms
from typing import List
from copy import deepcopy
from scipy.optimize import minimize_scalar
from scipy.spatial.distance import cdist


def ligand_metal_docking(
    ligand: Atoms,
    bonding_sites_clusters: List[List[int]],
    metal_center: Atoms,
    distance: float
) -> Atoms:
    """
    Places metal centers at each bonding site of the ligand, facing away from the ligand,
    such that the distance between each bonding site atom and the closest atom of the metal
    center is equal to the specified distance. The metal center's center of mass remains
    along the bonding site centroid vector.

    Parameters:
    - ligand: ASE Atoms object representing the ligand.
    - bonding_sites_clusters: List of clusters, each containing atom indices of a bonding site.
    - metal_center: ASE Atoms object representing the metal center to be placed.
    - distance: Distance in Angstroms between each bonding site atom and the closest metal center atom.

    Returns:
    - combined_atoms: ASE Atoms object with the metal centers added to the ligand.
    """

    # Copy the ligand
    ligand_metal_docked = ligand.copy()

    # Compute the center of mass of the ligand
    ligand_center_of_mass = ligand.get_center_of_mass()

    for bonding_site in bonding_sites_clusters:
        # Get positions of the bonding site atoms
        bonding_site_positions = ligand.positions[bonding_site]

        # Compute the centroid of the bonding site atoms
        bonding_site_centroid = np.mean(bonding_site_positions, axis=0)

        # Compute the vector from ligand's center of mass to bonding site centroid
        vector = bonding_site_centroid - ligand_center_of_mass

        # Normalize the vector
        unit_vector = vector / np.linalg.norm(vector)

        # Define the objective function
        def objective(scalar):
            # Compute the position along the centroid vector at distance scalar
            metal_center_position = ligand_center_of_mass + unit_vector * scalar

            # Translate the metal center so that its center of mass is at this position
            metal_atoms = metal_center.copy()
            metal_center_of_mass = metal_atoms.get_center_of_mass()
            translation = metal_center_position - metal_center_of_mass
            metal_atoms.translate(translation)

            # Compute distances from bonding site atoms to all atoms in the metal center
            distances = cdist(bonding_site_positions, metal_atoms.positions)

            # For each bonding site atom, find the minimum distance to the metal center atoms
            min_distances = np.min(distances, axis=1)

            # Compute the differences between these distances and the specified distance
            differences = min_distances - distance

            # Return the sum of squares of differences
            return np.sum(differences**2)

        # Use an optimizer to find the scalar that minimizes the objective function
        res = minimize_scalar(
            objective,
            bounds=(0, 10 * distance),  # Adjust upper bound as needed
            method='bounded'
        )

        optimal_scalar = res.x

        # Compute the final position of the metal center
        metal_center_position = ligand_center_of_mass + unit_vector * optimal_scalar

        # Translate the metal center to the final position
        metal_atoms = metal_center.copy()
        metal_center_of_mass = metal_atoms.get_center_of_mass()
        translation = metal_center_position - metal_center_of_mass
        metal_atoms.translate(translation)

        # Append the metal atoms to the combined_atoms
        ligand_metal_docked += metal_atoms

    return ligand_metal_docked
