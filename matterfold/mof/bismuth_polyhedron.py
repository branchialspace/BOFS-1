# Bismuth polyhedron metal center

from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
from ase.io import write
import numpy as np


def generate_fibonacci_sphere(samples: int, radius: float = 1.0) -> np.ndarray:
    """
    Generate points on a sphere using the Fibonacci lattice method.

    :param samples: Number of points to generate.
    :param radius: Radius of the sphere.
    :return: Array of shape (samples, 3) with 3D coordinates of points on the sphere.
    """
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius_at_y = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius_at_y
        z = np.sin(theta) * radius_at_y

        points.append((x * radius, y * radius, z * radius))

    return np.array(points)

def generate_bismuth_polyhedron(num_atoms: int) -> Atoms:
    """
    Generate a bismuth polyhedron with a specified number of atoms.

    :param num_atoms: Number of atoms.
    :return: ASE Atoms object with bismuth atoms arranged in a polyhedral structure.
    """
    if num_atoms < 2:
        raise ValueError("Number of atoms must be at least 2 to form a polyhedron.")

    # Generate positions on a sphere using the Fibonacci lattice method
    positions = generate_fibonacci_sphere(num_atoms)
    
    # Create the Atoms object
    atoms = Atoms(f'Bi{num_atoms}', positions=positions)
    
    # Set up the Lennard-Jones calculator
    lj_calc = LennardJones(sigma=3.5, epsilon=0.2)
    atoms.calc = lj_calc
    
    # Optimize the structure
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.01)

    filename = f'Bi{num_atoms}_polyhedron.xyz'
    write(filename, atoms)
    
    return atoms
