# Bismuth polyhedron metal center

import numpy as np
from scipy.spatial.distance import pdist, squareform
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
from ase.io import write


# Generate a close-packed cluster of spheres using a force-biased algorithm
def sphere_pack_cluster(n_atoms, radius=1.5, max_iterations=2000):
    # Initialize random positions in a cube
    positions = np.random.rand(n_atoms, 3) * (n_atoms**(1/3) * radius)
    
    for _ in range(max_iterations):
        # Calculate pairwise distances
        distances = squareform(pdist(positions))
        np.fill_diagonal(distances, np.inf)  # Avoid self-interactions
        
        # Calculate unit vectors between atoms
        diff = positions[:, np.newaxis] - positions
        with np.errstate(invalid='ignore', divide='ignore'):
            unit_vectors = diff / distances[:, :, np.newaxis]
        unit_vectors = np.nan_to_num(unit_vectors)
        
        # Calculate repulsive forces (inversely proportional to distance)
        forces = np.maximum(0, 2*radius - distances)[:, :, np.newaxis] * unit_vectors
        
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

# Generate a bismuth polyhedron with a specified number of atoms
def generate_bismuth_polyhedron(num_atoms: int) -> Atoms:
    if num_atoms == 1:
        positions = np.array([[0.0, 0.0, 0.0]])
    elif num_atoms >= 2:
        # Generate packed positions using the sphere packing method
        packed_positions = sphere_pack_cluster(num_atoms, radius=1.5)
        
        # Normalize to fit within a unit sphere
        max_distance = np.max(np.linalg.norm(packed_positions, axis=1))
        packed_positions /= max_distance
        
        # Scale the positions to a reasonable size for bismuth atoms (e.g., 3 Angstrom radius)
        positions = packed_positions * 3.0
    else:
        raise ValueError("Number of atoms must be at least 1")
    
    # Create the Atoms object
    atoms = Atoms(f'Bi{num_atoms}', positions=positions)
    
    # Set up the Lennard-Jones calculator
    lj_calc = LennardJones(sigma=3.2, epsilon=0.2)
    atoms.calc = lj_calc
    
    # Optimize the structure
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.01)
    
    filename = f'Bi{num_atoms}_polyhedron.xyz'
    write(filename, atoms)
    
    return atoms
