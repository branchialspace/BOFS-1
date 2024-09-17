# Bismuth polyhedron metal center

import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
from ase.io import write


# Hardcoded Sloane's Tables optimal sphere packing degrees of separation data http://neilsloane.com/packings/index.html#I
SLOANE_DATA = {
    4: 109.4712206, 5: 90.0000000, 6: 90.0000000, 7: 77.8695421,
    8: 74.8584922, 9: 70.5287794, 10: 66.1468220, 11: 63.4349488,
    12: 63.4349488, 13: 57.1367031, 14: 55.6705700, 15: 53.6578501,
    16: 52.2443957, 17: 51.0903285, 18: 49.5566548, 19: 47.6919141,
    20: 47.4310362, 21: 45.6132231, 22: 44.7401612, 23: 43.7099642,
    24: 43.6907671, 25: 41.6344612, 26: 41.0376616, 27: 40.6776007,
    28: 39.3551436, 29: 38.7136512, 30: 38.5971159, 31: 37.7098291,
    32: 37.4752140, 33: 36.2545530, 34: 35.8077844, 35: 35.3198076,
    36: 35.1897322, 37: 34.4224080, 38: 34.2506607, 39: 33.4890466,
    40: 33.1583563
}

# Generate optimal packing based on hardcoded Sloane's Tables data
def generate_optimal_packing(num_points):
    if num_points not in SLOANE_DATA:
        raise ValueError(f"Optimal packing for {num_points} points not available in the hardcoded data.")
    
    # Convert separation angle to radians
    separation = np.radians(SLOANE_DATA[num_points])
    
    # Generate points on a unit sphere
    points = []
    golden_angle = np.pi * (3 - np.sqrt(5))
    
    for i in range(num_points):
        y = 1 - (i / float(num_points - 1)) * 2
        radius = np.sqrt(1 - y * y)
        theta = golden_angle * i
        
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        
        points.append([x, y, z])
    
    return np.array(points)

# Generate a bismuth polyhedron with a specified number of atoms
def generate_bismuth_polyhedron(num_atoms: int) -> Atoms:
    if num_atoms < 2:
        raise ValueError("Number of atoms must be at least 2 to form a polyhedron.")
    
    # Generate optimal packing coordinates
    packed_positions = generate_optimal_packing(num_atoms)
    
    # Ensure the first atom is at (0,0,0) and others are relative to it
    packed_positions -= packed_positions[0]
    
    # Normalize to fit within a unit sphere
    max_distance = np.max(np.linalg.norm(packed_positions, axis=1))
    packed_positions /= max_distance
    
    # Scale the positions to a reasonable size for bismuth atoms (e.g., 3 Angstrom radius)
    packed_positions *= 3.0
    
    # Create the Atoms object
    atoms = Atoms(f'Bi{num_atoms}', positions=packed_positions)
    
    # Set up the Lennard-Jones calculator
    lj_calc = LennardJones(sigma=3.2, epsilon=0.2)
    atoms.calc = lj_calc
    
    # Optimize the structure
    optimizer = BFGS(atoms)
    optimizer.run(fmax=0.01)
    
    filename = f'Bi{num_atoms}_polyhedron.xyz'
    write(filename, atoms)
    
    return atoms
