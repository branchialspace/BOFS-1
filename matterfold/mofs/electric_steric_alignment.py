# Bisected bonding-site charge symmetry for 1d bridging ligands
# Steric hindrance in single-cluster octahedral coordination
    # Bonding site angle = max distance from internal ligand bonds

import numpy as np
from scipy.optimize import minimize
from ase import Atoms

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
    weights = np.array([bonding_sites[i + 1]['non_bonded_electrons'] for i in range(num_atoms)])

    # Normalize weights to avoid division by zero
    total_weight = np.sum(weights)
    if total_weight == 0:
        raise ValueError("Total unbonded electrons in the ligand are zero.")

    # Define the objective function to maximize the distance between charge centers
    def objective_function(params):
        theta, phi = params
        n = np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ])

        # Project positions onto the normal vector
        projections = positions @ n

        # Sort projections and corresponding weights
        sorted_indices = np.argsort(projections)
        sorted_projections = projections[sorted_indices]
        sorted_weights = weights[sorted_indices]

        # Compute cumulative weights
        cumulative_weights = np.cumsum(sorted_weights)
        total_weight = cumulative_weights[-1]

        # Find the plane position that balances the weights
        idx = np.searchsorted(cumulative_weights, total_weight / 2)
        d_candidates = sorted_projections[idx - 1:idx + 1] if idx > 0 else sorted_projections[:1]
        best_D = -np.inf
        best_d = None

        for d in d_candidates:
            # Assign atoms to each side of the plane
            side = np.where(projections >= d, 1, -1)
            W_plus = np.sum(weights[side == 1])
            W_minus = np.sum(weights[side == -1])

            # Calculate centers of unbonded electron densities
            if W_plus > 0 and W_minus > 0:
                C_plus = np.sum(positions[side == 1] * weights[side == 1, np.newaxis], axis=0) / W_plus
                C_minus = np.sum(positions[side == -1] * weights[side == -1, np.newaxis], axis=0) / W_minus

                # Compute distance between centers
                D = np.linalg.norm(C_plus - C_minus)

                # Check if this division is better
                weight_diff = abs(W_plus - W_minus)
                if weight_diff <= total_weight * 0.1 and D > best_D:  # Allowing 10% imbalance
                    best_D = D
                    best_d = d

        # Return negative distance to convert maximization to minimization
        return -best_D if best_D != -np.inf else np.inf

    # Initial guess for theta and phi
    initial_guess = [np.pi / 2, 0]

    # Bounds for theta and phi
    bounds = [(0, np.pi), (0, 2 * np.pi)]

    # Optimize using scipy's minimize function
    result = minimize(
        objective_function,
        initial_guess,
        bounds=bounds,
        method='L-BFGS-B',
        options={'ftol': 1e-6, 'disp': True}
    )

    # Extract the best normal vector
    theta_opt, phi_opt = result.x
    n_opt = np.array([
        np.sin(theta_opt) * np.cos(phi_opt),
        np.sin(theta_opt) * np.sin(phi_opt),
        np.cos(theta_opt)
    ])

    # Find the best plane position d
    projections = positions @ n_opt
    sorted_indices = np.argsort(projections)
    sorted_projections = projections[sorted_indices]
    sorted_weights = weights[sorted_indices]
    cumulative_weights = np.cumsum(sorted_weights)
    total_weight = cumulative_weights[-1]
    idx = np.searchsorted(cumulative_weights, total_weight / 2)
    d_candidates = sorted_projections[idx - 1:idx + 1] if idx > 0 else sorted_projections[:1]
    best_D = -np.inf
    best_d = None

    for d in d_candidates:
        side = np.where(projections >= d, 1, -1)
        W_plus = np.sum(weights[side == 1])
        W_minus = np.sum(weights[side == -1])

        if W_plus > 0 and W_minus > 0:
            C_plus = np.sum(positions[side == 1] * weights[side == 1, np.newaxis], axis=0) / W_plus
            C_minus = np.sum(positions[side == -1] * weights[side == -1, np.newaxis], axis=0) / W_minus
            D = np.linalg.norm(C_plus - C_minus)
            weight_diff = abs(W_plus - W_minus)
            if weight_diff <= total_weight * 0.1 and D > best_D:
                best_D = D
                best_d = d

    if best_d is None:
        raise ValueError("Could not find a suitable plane to bisect the ligand.")

    # Return the plane parameters
    best_plane = {
        'normal_vector': n_opt,
        'point_on_plane': n_opt * best_d
    }

    return best_plane
