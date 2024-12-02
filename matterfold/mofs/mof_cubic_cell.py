# MOF primitive cubic lattice unit cell

import numpy as np
from ase.io import write
from ase import Atoms
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation as R


def mof_cell(
        combined_structure: Atoms,
        metal_center: Atoms,
        ligand: Atoms,
        bonding_sites: list,
    ) -> Atoms:
    """
    Creates a primitive cubic unit cell of our MOF model.
    Uses relative positions from the combined_structure to coordinate ligands to surface atoms of the cluster.
    
    Parameters:
    - combined_structure: ASE Atoms object of a linear bridging ligand with metal centers placed at each bonding site
    - metal_center: ASE Atoms object of the metal cluster
    - ligand: ASE Atoms object of the ligand
    - bonding_sites: List of lists, each containing atom indices (0-based) of a bonding site

    Returns:
    - unit_cell: ASE Atoms object of our MOF primitive cubic lattice unit cell

    """
    # Number of metal atoms and ligand atoms
    n_metal_single = len(metal_center)  # Should be 6
    n_ligand = len(ligand)
    
    # Extract ligand and metal positions from combined_structure
    ligand_positions_combined = combined_structure.positions[:n_ligand]
    metal_positions_combined = combined_structure.positions[-2*n_metal_single:]  # Get all metal positions
    
    # Split into two metal centers (each should have 6 atoms)
    metal_center1 = metal_positions_combined[:n_metal_single]
    metal_center2 = metal_positions_combined[n_metal_single:]
    
    # Calculate metal center centroids
    centroid1 = np.mean(metal_center1, axis=0)
    centroid2 = np.mean(metal_center2, axis=0)
    
    # Find the first bonding site centroid
    bonding_site_positions = [ligand_positions_combined[idx] for idx in bonding_sites[0]]
    bonding_site_centroid = np.mean(bonding_site_positions, axis=0)
    
    # Determine which metal center is closer to the bonding site
    dist1 = np.linalg.norm(bonding_site_centroid - centroid1)
    dist2 = np.linalg.norm(bonding_site_centroid - centroid2)
    
    # Select the closer metal center and its positions
    metal_positions = metal_center1 if dist1 < dist2 else metal_center2
    metal_centroid = centroid1 if dist1 < dist2 else centroid2
        
    # Find coordinating surface atom in the selected metal center
    metal_hull = ConvexHull(metal_positions)
    surface_atom_indices = np.unique(metal_hull.simplices.flatten())
    
    # Find the surface atom closest to the bonding site
    min_distance = float('inf')
    surface_atom_idx = None
    
    for idx in surface_atom_indices:
        surface_pos = metal_positions[idx]
        distance = np.linalg.norm(bonding_site_centroid - surface_pos)
        if distance < min_distance:
            min_distance = distance
            surface_atom_idx = idx
            
    surface_atom_pos = metal_positions[surface_atom_idx]
    
    # Calculate relative positions of ligand atoms to the surface atom
    ligand_relative_positions = ligand_positions_combined - surface_atom_pos
    
    # Calculate the original ligand direction (mean direction)
    original_direction = ligand_relative_positions.mean(axis=0)
    original_direction /= np.linalg.norm(original_direction)
    
    # Calculate maximum ligand extent from coordination point
    max_ligand_extent = np.max(np.linalg.norm(ligand_relative_positions, axis=1))
    
    # Calculate metal center extent in x, y, z directions from its centroid
    metal_positions_centered = metal_center.positions - np.mean(metal_center.positions, axis=0)
    metal_extent = np.max(np.abs(metal_positions_centered))
    
    # Calculate lattice constant based on metal center extent plus ligand extent in each direction
    lattice_constant = 2 * (metal_extent + max_ligand_extent)
    
    # Create unit cell with centered metal center
    unit_cell = Atoms(
        symbols=metal_center.get_chemical_symbols(),
        positions=metal_center.get_positions(),
        cell=[lattice_constant, lattice_constant, lattice_constant],
        pbc=True
    )
    unit_cell.center()
    
    # Get positions after centering
    unit_cell_positions = unit_cell.get_positions()
    center_pos = unit_cell.get_center_of_mass()
    
    # Define directions for x, y, z
    directions = [
        np.array([1, 0, 0]),  # +x-direction
        np.array([0, 1, 0]),  # +y-direction
        np.array([0, 0, 1])   # +z-direction
    ]
    
    surface_atoms = []
    for direction in directions:
        # Project positions onto direction vector
        projections = np.dot(unit_cell_positions - center_pos, direction)
        surface_idx = np.argmax(projections)
        surface_atoms.append(surface_idx)
    
    # Add coordinated ligands to each surface atom with correct orientation
    for surface_idx in surface_atoms:
        new_surface_pos = unit_cell_positions[surface_idx]
        
        # Calculate the new direction vector from metal centroid to surface atom
        new_direction = new_surface_pos - center_pos
        new_direction /= np.linalg.norm(new_direction)  # Normalize
        
        # Calculate rotation needed to align original_direction to new_direction
        rotation_axis = np.cross(original_direction, new_direction)
        norm_axis = np.linalg.norm(rotation_axis)
        
        if norm_axis < 1e-6:
            # Directions are the same or opposite
            if np.dot(original_direction, new_direction) > 0:
                rotation_matrix = np.identity(3)  # No rotation needed
            else:
                # 180-degree rotation around an orthogonal axis
                if not np.allclose(original_direction, [1, 0, 0]):
                    orthogonal = np.array([1, 0, 0])
                else:
                    orthogonal = np.array([0, 1, 0])
                rotation_axis = np.cross(original_direction, orthogonal)
                rotation_axis /= np.linalg.norm(rotation_axis)
                rotation = R.from_rotvec(np.pi * rotation_axis)
                rotation_matrix = rotation.as_matrix()
        else:
            # Calculate the angle between the original and new directions
            rotation_axis /= norm_axis
            angle = np.arccos(np.clip(np.dot(original_direction, new_direction), -1.0, 1.0))
            rotation = R.from_rotvec(rotation_axis * angle)
            rotation_matrix = rotation.as_matrix()
        
        # Apply rotation to ligand_relative_positions
        rotated_ligand_rel_pos = np.dot(ligand_relative_positions, rotation_matrix.T)
        
        # Create new ligand with rotated positions relative to new surface atom
        ligand_copy = ligand.copy()
        new_ligand_positions = rotated_ligand_rel_pos + new_surface_pos
        ligand_copy.set_positions(new_ligand_positions)
        
        unit_cell += ligand_copy
    
    metal_formula = metal_center.get_chemical_formula()
    ligand_formula = ligand.get_chemical_formula()
    filename = f"{metal_formula}_{ligand_formula}_cubic_cell.xyz"
    write(filename, unit_cell)
    
    return unit_cell
