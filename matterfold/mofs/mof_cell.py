# MOF primitive cubic/ rhombohedral lattice unit cell

from ase import Atoms, Atom
from ase.io import write
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

def build_metal_organic_crystal(
    combined_structure: Atoms,
    ligand: Atoms,
    metal_center: Atoms,
    bonding_sites: list,
    tolerance: float = 0.1
) -> Atoms:
    """
    Constructs a primitive cubic/rhombohedral crystal structure from metal-ligand building blocks.
    
    Parameters:
    - combined_structure: ASE Atoms object from ligand_metal_docking containing single M-L-M unit
    - ligand: Original ligand ASE Atoms object
    - metal_center: Original metal cluster ASE Atoms object
    - bonding_sites: List of lists containing atom indices of bonding sites
    - tolerance: Tolerance for geometric comparisons (default: 0.1 Å)
    
    Returns:
    - crystal: ASE Atoms object with complete unit cell and periodic boundary conditions
    """
    # 1. Structure Analysis
    positions = combined_structure.get_positions()
    symbols = combined_structure.get_chemical_symbols()
    
    # Identify metal cluster positions by matching metal symbols
    metal_symbol = metal_center[0].symbol
    metal_indices = [i for i, sym in enumerate(symbols) if sym == metal_symbol]
    
    # Group metal atoms into clusters
    n_metal_atoms = len(metal_center)
    metal_clusters = [metal_indices[i:i+n_metal_atoms] 
                     for i in range(0, len(metal_indices), n_metal_atoms)]
    
    # Extract reference metal cluster information
    reference_clusters = []
    for cluster_indices in metal_clusters:
        cluster_pos = positions[cluster_indices]
        centroid = np.mean(cluster_pos, axis=0)
        
        # Find coordinating atom using ConvexHull
        hull = ConvexHull(cluster_pos)
        coord_atom_idx = cluster_indices[np.unique(hull.simplices.flatten())[0]]
        coord_vector = positions[coord_atom_idx] - centroid
        coord_vector = coord_vector / np.linalg.norm(coord_vector)
        
        # Calculate full orientation matrix
        ref_pos = metal_center.get_positions()
        ref_centroid = np.mean(ref_pos, axis=0)
        
        reference_clusters.append({
            'centroid': centroid,
            'positions': cluster_pos,
            'coord_vector': coord_vector,
            'relative_positions': cluster_pos - centroid
        })
    
    # Calculate cell vector length (metal-to-metal distance)
    cell_length = np.linalg.norm(reference_clusters[1]['centroid'] - 
                                reference_clusters[0]['centroid'])
    
    # Calculate angle between coordination vectors
    coord_angle = np.arccos(np.dot(reference_clusters[0]['coord_vector'],
                                  reference_clusters[1]['coord_vector']))
    is_rhombohedral = abs(coord_angle - np.pi) > tolerance
    
    # 3. Unit Cell Construction
    if is_rhombohedral:
        # Calculate rhombohedral angles
        cos_alpha = np.dot(reference_clusters[0]['coord_vector'], 
                          reference_clusters[1]['coord_vector'])
        alpha = np.arccos(cos_alpha)
        cell = np.array([
            [cell_length, 0, 0],
            [cell_length * np.cos(alpha), cell_length * np.sin(alpha), 0],
            [cell_length * np.cos(alpha), 
             cell_length * np.cos(alpha) * np.cos(alpha),
             cell_length * np.sqrt(1 - 2*np.cos(alpha)**2 + np.cos(alpha)**2)]
        ])
    else:
        # Cubic cell
        cell = np.eye(3) * cell_length
    
    # 4. Create Crystal Structure
    crystal = Atoms()
    
    def get_cluster_orientation(position, direction_vectors):
        """Calculate cluster orientation for given position and connection directions"""
        # Use reference cluster closest to desired orientation
        ref_idx = 0 if len(direction_vectors) == 1 else (
            0 if np.dot(reference_clusters[0]['coord_vector'], direction_vectors[0]) > 
               np.dot(reference_clusters[1]['coord_vector'], direction_vectors[0]) else 1
        )
        
        ref_cluster = reference_clusters[ref_idx]
        
        # Calculate rotation to align coordination vector with first direction
        v1 = ref_cluster['coord_vector']
        v2 = direction_vectors[0]
        
        # Calculate rotation axis and angle
        rotation_axis = np.cross(v1, v2)
        if np.any(rotation_axis):
            rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
            angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
            
            # Create rotation matrix
            c = np.cos(angle)
            s = np.sin(angle)
            x, y, z = rotation_axis
            R = np.array([[c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
                         [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
                         [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])
            
            rotated_positions = ref_cluster['relative_positions'] @ R.T
        else:
            rotated_positions = ref_cluster['relative_positions']
            
        return rotated_positions + position
    
    # Add metal clusters at lattice points
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pos = np.dot([i, j, k], cell)
                
                # Determine connection directions for this cluster
                directions = []
                if i < 1:
                    directions.append(cell[0] / np.linalg.norm(cell[0]))
                if i > 0:
                    directions.append(-cell[0] / np.linalg.norm(cell[0]))
                if j < 1:
                    directions.append(cell[1] / np.linalg.norm(cell[1]))
                if j > 0:
                    directions.append(-cell[1] / np.linalg.norm(cell[1]))
                if k < 1:
                    directions.append(cell[2] / np.linalg.norm(cell[2]))
                if k > 0:
                    directions.append(-cell[2] / np.linalg.norm(cell[2]))
                
                # Calculate cluster positions with proper orientation
                cluster_positions = get_cluster_orientation(pos, directions)
                
                # Add oriented cluster
                cluster_atoms = Atoms([metal_center[0].symbol] * len(metal_center),
                                    positions=cluster_positions)
                crystal += cluster_atoms
    
    # Add connecting ligands
    for i in range(2):
        for j in range(2):
            for k in range(2):
                if i < 1 or j < 1 or k < 1:  # Avoid double-counting edges
                    pos1 = np.dot([i, j, k], cell)
                    
                    if i < 1:
                        pos2 = np.dot([i+1, j, k], cell)
                        add_ligand(crystal, ligand, pos1, pos2, cell[0], bonding_sites)
                    if j < 1:
                        pos2 = np.dot([i, j+1, k], cell)
                        add_ligand(crystal, ligand, pos1, pos2, cell[1], bonding_sites)
                    if k < 1:
                        pos2 = np.dot([i, j, k+1], cell)
                        add_ligand(crystal, ligand, pos1, pos2, cell[2], bonding_sites)
    
    # 5. Set Periodic Boundary Conditions
    crystal.set_cell(cell)
    crystal.set_pbc(True)
    
    # Write structure to file
    write('metal_organic_crystal.xyz', crystal)
    
    return crystal

def add_ligand(crystal, ligand, pos1, pos2, direction, bonding_sites):
    """Helper function to add a ligand between two metal centers"""
    new_ligand = ligand.copy()
    
    # Calculate translation and rotation
    mid_point = (pos1 + pos2) / 2
    direction = direction / np.linalg.norm(direction)
    
    # Align ligand with direction
    current_axis = np.array([1, 0, 0])  # Assume ligand aligned along x-axis initially
    rotation_axis = np.cross(current_axis, direction)
    
    if np.any(rotation_axis):
        angle = np.arccos(np.dot(current_axis, direction))
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        new_ligand.rotate(angle * 180 / np.pi, rotation_axis)
    
    new_ligand.translate(mid_point)
    crystal += new_ligand
