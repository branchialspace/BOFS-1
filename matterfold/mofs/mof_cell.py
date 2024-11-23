# MOF primitive cubic/ rhombohedral lattice unit cell

from ase import Atoms
from ase.io import write
import numpy as np


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
    
    # Extract reference metal clusters and connecting vectors
    reference_clusters = []
    for cluster_indices in metal_clusters:
        cluster_pos = positions[cluster_indices]
        centroid = np.mean(cluster_pos, axis=0)
        
        # Calculate reference connecting vector using bonding sites
        site_positions = [positions[i] for i in bonding_sites[0]]
        connecting_point = np.mean(site_positions, axis=0)
        connect_vector = connecting_point - centroid
        connect_vector = connect_vector / np.linalg.norm(connect_vector)
        
        reference_clusters.append({
            'centroid': centroid,
            'positions': cluster_pos,
            'relative_positions': cluster_pos - centroid,
            'connect_vector': connect_vector
        })
    
    # Calculate cell parameters
    cell_vector = reference_clusters[1]['centroid'] - reference_clusters[0]['centroid']
    cell_length = np.linalg.norm(cell_vector)
    
    # Define cubic cell
    cell = np.eye(3) * cell_length
    
    # 4. Create Crystal Structure
    crystal = Atoms()
    
    def get_cluster_orientation(position, directions):
        """Calculate cluster orientation for given position and connection directions"""
        if not directions:
            return reference_clusters[0]['relative_positions'] + position
            
        # Get reference orientation from combined structure
        ref_connect = reference_clusters[0]['connect_vector']
        
        # Create orthogonal basis from first direction
        v1 = directions[0]
        v1 = v1 / np.linalg.norm(v1)
        
        # Find rotation that aligns reference vector with first direction
        rotation_axis = np.cross(ref_connect, v1)
        if np.any(rotation_axis):
            rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
            angle = np.arccos(np.clip(np.dot(ref_connect, v1), -1.0, 1.0))
            
            # Create rotation matrix
            c = np.cos(angle)
            s = np.sin(angle)
            x, y, z = rotation_axis
            R = np.array([[c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
                         [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
                         [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])
        else:
            R = np.eye(3)
            
        # Apply rotation to relative positions
        rotated_positions = np.dot(reference_clusters[0]['relative_positions'], R.T)
        
        return rotated_positions + position
    
    # Add metal clusters at lattice points
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pos = np.dot([i, j, k], cell)
                
                # Determine connection directions
                directions = []
                if i < 1:
                    directions.append(cell[0] / cell_length)
                if i > 0:
                    directions.append(-cell[0] / cell_length)
                if j < 1:
                    directions.append(cell[1] / cell_length)
                if j > 0:
                    directions.append(-cell[1] / cell_length)
                if k < 1:
                    directions.append(cell[2] / cell_length)
                if k > 0:
                    directions.append(-cell[2] / cell_length)
                
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
                        pos2 = pos1 + cell[0]
                        add_ligand_with_orientation(crystal, ligand, pos1, pos2, 
                                                 reference_clusters[0], bonding_sites)
                    if j < 1:
                        pos2 = pos1 + cell[1]
                        add_ligand_with_orientation(crystal, ligand, pos1, pos2, 
                                                 reference_clusters[0], bonding_sites)
                    if k < 1:
                        pos2 = pos1 + cell[2]
                        add_ligand_with_orientation(crystal, ligand, pos1, pos2, 
                                                 reference_clusters[0], bonding_sites)
    
    # 5. Set Periodic Boundary Conditions
    crystal.set_cell(cell)
    crystal.set_pbc(True)

    # Write structure to file
    write('metal_organic_crystal.xyz', crystal)
    
    return crystal

def add_ligand_with_orientation(crystal, ligand, pos1, pos2, ref_cluster, bonding_sites):
    """Add ligand between two points with orientation matching the reference structure"""
    new_ligand = ligand.copy()
    
    # Get reference connecting vector and current direction
    ref_connect = ref_cluster['connect_vector']
    current_direction = pos2 - pos1
    current_direction = current_direction / np.linalg.norm(current_direction)
    
    # Calculate rotation to align with current direction
    rotation_axis = np.cross(ref_connect, current_direction)
    if np.any(rotation_axis):
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        angle = np.arccos(np.clip(np.dot(ref_connect, current_direction), -1.0, 1.0))
        
        # Create rotation matrix
        c = np.cos(angle)
        s = np.sin(angle)
        x, y, z = rotation_axis
        R = np.array([[c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
                     [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
                     [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])
        
        # Apply rotation
        positions = new_ligand.get_positions()
        rotated_positions = np.dot(positions, R.T)
        new_ligand.set_positions(rotated_positions)
    
    # Place ligand at midpoint
    mid_point = (pos1 + pos2) / 2
    new_ligand.translate(mid_point)
    
    crystal += new_ligand
