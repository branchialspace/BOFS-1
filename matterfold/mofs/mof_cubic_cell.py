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
    """
    # Number of metal atoms and ligand atoms
    n_metal = len(metal_center)
    n_ligand = len(ligand)
    
    # Extract ligand and metal positions from combined_structure
    ligand_positions_combined = combined_structure.positions[:n_ligand]
    metal_positions_combined = combined_structure.positions[-n_metal:]
    metal_centroid = np.mean(metal_positions_combined, axis=0)
    
    # Find coordinating surface atom in combined_structure
    metal_hull = ConvexHull(metal_positions_combined)
    surface_atom_indices_combined = np.unique(metal_hull.simplices.flatten())
    
    # Assume that the first surface atom is the coordinating one for ligand
    surface_atom_idx_combined = surface_atom_indices_combined[0]
    surface_atom_pos_combined = metal_positions_combined[surface_atom_idx_combined]
    
    # Calculate relative positions of ligand atoms to the surface atom in combined_structure
    ligand_relative_positions = ligand_positions_combined - surface_atom_pos_combined
    
    # Calculate the original ligand direction (mean direction)
    original_direction = ligand_relative_positions.mean(axis=0)
    original_direction /= np.linalg.norm(original_direction)
    
    # Calculate lattice constant based on the combined structure dimensions
    max_distance = np.max(np.linalg.norm(ligand_relative_positions, axis=1))
    lattice_constant = max_distance * 6  # Ensure enough space
    
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
    filename = f"{metal_formula}_{ligand_formula}_cell.xyz"
    write(filename, unit_cell)
    
    return unit_cell
