# Metal polyhedron cluster as metal center building unit

import numpy as np
from scipy.spatial.distance import pdist, squareform
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.io import write


# Generate a close-packed cluster of spheres using a force-biased algorithm
def sphere_pack_cluster(species, n_atoms, max_iterations=2000, seed=42):
    # Get covalent radius of species
    radius = covalent_radii[species]
    # Initialize random positions in a cube, random seed for reproducibility
    np.random.seed(seed)
    positions = np.random.rand(n_atoms, 3) * (n_atoms**(1/3) * radius)

    for _ in range(max_iterations):
        # Calculate pairwise distances
        distances = squareform(pdist(positions))
        np.fill_diagonal(distances, np.inf)  # Avoid self-interactions

        # Calculate unit vectors between atoms
        diff = positions[:, np.newaxis] - positions
        unit_vectors = diff / distances[:, :, np.newaxis]
        unit_vectors = np.nan_to_num(unit_vectors)

        # Calculate repulsive forces (inversely proportional to distance)
        # radius multiplied by 2.2 rather than 2 to better estimate metallic radii
        forces = np.maximum(0, 2.2*radius - distances)[:, :, np.newaxis] * unit_vectors

        # Sum forces on each atom
        total_forces = np.sum(forces, axis=1)

        # Move atoms based on forces
        positions += total_forces * 0.1  # Small step size for stability

        # Apply weak central force to keep cluster together
        center = np.mean(positions, axis=0)
        to_center = center - positions
        positions += to_center * 0.01

    # Center the cluster at the origin
    positions -= np.mean(positions, axis=0)

    return positions

# Generate a metal polyhedron with a specified number of atoms
def generate_metal_polyhedron(species: str, num_atoms: int) -> Atoms:
    if num_atoms == 1:
        positions = np.array([[0.0, 0.0, 0.0]])
    elif num_atoms >= 2:
        # Generate packed positions using the sphere packing method
        species_str = species.capitalize()
        species = atomic_numbers[species_str]
        packed_positions = sphere_pack_cluster(species, num_atoms)
        positions = packed_positions
    else:
        raise ValueError("Number of atoms must be at least 1")

    # Create the Atoms object
    atoms = Atoms(f'{species_str}{num_atoms}', positions=positions)

    filename = f'{species_str}{num_atoms}_polyhedron.xyz'
    write(filename, atoms)

    return atoms
