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
    n_metal_single = len(metal_center)
    n_ligand = len(ligand)
    
    # Extract ligand and metal positions from combined_structure
    ligand_positions_combined = combined_structure.positions[:n_ligand]
    metal_positions_combined = combined_structure.positions[-2*n_metal_single:]  # Get all metal positions
    
    # Split into two metal centers
    metal_center1 = metal_positions_combined[:n_metal_single]
    metal_center2 = metal_positions_combined[n_metal_single:]
    
    # Calculate metal center centroids
    centroid1 = np.mean(metal_center1, axis=0)
    centroid2 = np.mean(metal_center2, axis=0)
    
    # Find the first bonding site centroid
    bonding_site_positions = [ligand_positions_combined[idx] for idx in bonding_sites[0]]
    bonding_site_centroid = np.mean(bonding_site_positions, axis=0)
    
    # Find coordinating metal atoms from both metal centers
    surface_atom1_idx = None
    surface_atom2_idx = None
    min_distance1 = float('inf')
    min_distance2 = float('inf')
    
    # Find surface atoms for first metal center
    hull1 = ConvexHull(metal_center1)
    surface_indices1 = np.unique(hull1.simplices.flatten())
    for idx in surface_indices1:
        surface_pos = metal_center1[idx]
        distance = np.linalg.norm(bonding_site_centroid - surface_pos)
        if distance < min_distance1:
            min_distance1 = distance
            surface_atom1_idx = idx
            
    # Find surface atoms for second metal center
    hull2 = ConvexHull(metal_center2)
    surface_indices2 = np.unique(hull2.simplices.flatten())
    for idx in surface_indices2:
        surface_pos = metal_center2[idx]
        distance = np.linalg.norm(bonding_site_centroid - surface_pos)
        if distance < min_distance2:
            min_distance2 = distance
            surface_atom2_idx = idx
    
    # Get positions of coordinating metal atoms
    coord_metal1_pos = metal_center1[surface_atom1_idx]
    coord_metal2_pos = metal_center2[surface_atom2_idx]
    
    # Calculate distance between coordinating metal atoms
    metal_metal_distance = np.linalg.norm(coord_metal2_pos - coord_metal1_pos)
    
    # Determine which metal center is closer to the bonding site
    dist1 = np.linalg.norm(bonding_site_centroid - centroid1)
    dist2 = np.linalg.norm(bonding_site_centroid - centroid2)
    
    # Select the closer metal center and its positions
    if dist1 < dist2:
        metal_positions = metal_center1
        metal_centroid = centroid1
        surface_atom_pos = coord_metal1_pos
        unused_metal_center = metal_center2
        unused_coord_metal_pos = coord_metal2_pos
    else:
        metal_positions = metal_center2
        metal_centroid = centroid2
        surface_atom_pos = coord_metal2_pos
        unused_metal_center = metal_center1
        unused_coord_metal_pos = coord_metal1_pos
        
    # Calculate relative positions of ligand atoms to the surface atom
    ligand_relative_positions = ligand_positions_combined - surface_atom_pos
    
    # Calculate the original ligand direction (mean direction)
    original_direction = ligand_relative_positions.mean(axis=0)
    original_direction /= np.linalg.norm(original_direction)
    
    # Find surface atoms for the metal center we are using
    hull = ConvexHull(metal_positions)
    surface_indices = np.unique(hull.simplices.flatten())
    max_distance = -float('inf')
    opposite_surface_atom_idx = None
    for idx in surface_indices:
        pos = metal_positions[idx]
        distance = np.linalg.norm(pos - surface_atom_pos)
        if distance > max_distance:
            max_distance = distance
            opposite_surface_atom_idx = idx
    opposite_surface_atom_pos = metal_positions[opposite_surface_atom_idx]
    
    # Calculate lattice constant based on the distance between the opposite surface atom
    # and the coordinating atom of the unused metal center
    lattice_constant = np.linalg.norm(opposite_surface_atom_pos - unused_coord_metal_pos)

    
    # Create unit cell
    unit_cell = Atoms(
        symbols=metal_center.get_chemical_symbols(),
        positions=metal_center.get_positions(),
        cell=[lattice_constant, lattice_constant, lattice_constant],
        pbc=True
    )
    
    # Get positions
    unit_cell_positions = unit_cell.get_positions()
    center_pos = unit_cell.get_center_of_mass()
    
    # Define directions for x, y, z
    directions = [
        np.array([1, 0, 0]),  # +x-direction
        np.array([0, 1, 0]),  # +y-direction
        np.array([0, 0, 1])   # +z-direction
    ]
    
    surface_atoms = []
    ligand_rotations = []  # To store rotation matrices for each ligand
    ligand_translations = []  # To store translations (surface positions) for each ligand
    far_side_coord_atom_positions = []  # Store far-side coordinating atom positions
    vectors_v = []  # Vectors from surface atom to far-side coordinating atom
    unit_vectors_u = []  # Target vectors along lattice vectors
    
    # Compute vector from bonding site centroid to unused coord_metal_pos
    unused_coord_vector = unused_coord_metal_pos - bonding_site_centroid
    
    for idx, direction in enumerate(directions):
        # Project positions onto direction vector
        projections = np.dot(unit_cell_positions - center_pos, direction)
        surface_idx = np.argmax(projections)
        surface_atoms.append(surface_idx)
        
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
        
        # Store rotation and translation
        ligand_rotations.append(rotation_matrix)
        ligand_translations.append(new_surface_pos)
        
        # Apply rotation to ligand_relative_positions
        rotated_ligand_rel_pos = np.dot(ligand_relative_positions, rotation_matrix.T)
        
        # Create new ligand with rotated positions relative to new surface atom
        ligand_copy = ligand.copy()
        new_ligand_positions = rotated_ligand_rel_pos + new_surface_pos
        ligand_copy.set_positions(new_ligand_positions)
        
        unit_cell += ligand_copy
        
        # Compute far-side coordinating atom position
        far_side_coord_atom_pos = np.dot(unused_coord_vector, rotation_matrix.T) + new_surface_pos
        far_side_coord_atom_positions.append(far_side_coord_atom_pos)
        
        # Compute vector v_i from surface atom to far-side coordinating atom
        v_i = far_side_coord_atom_pos - new_surface_pos
        vectors_v.append(v_i)
        
        # Target vector along lattice vector direction
        unit_vectors_u.append(direction * lattice_constant)
    
    # Adjust the position and rotation of the entire structure to align with unit cell boundaries
    
    # Stack vectors_v and unit_vectors_u
    V = np.array(vectors_v).T  # Shape (3, 3)
    U = np.array(unit_vectors_u).T  # Shape (3, 3)
    
    # Compute the covariance matrix
    H = V @ U.T
    
    # Compute SVD
    U_SVD, S, Vt_SVD = np.linalg.svd(H)
    
    # Compute rotation matrix
    R_opt = Vt_SVD.T @ U_SVD.T
    
    # Correct for reflection
    if np.linalg.det(R_opt) < 0:
        Vt_SVD[-1, :] *= -1
        R_opt = Vt_SVD.T @ U_SVD.T
    
    # Apply rotation to the entire unit cell
    positions = unit_cell.get_positions() - center_pos  # Center positions
    rotated_positions = np.dot(positions, R_opt.T)
    unit_cell.set_positions(rotated_positions + center_pos)
    
    # Update far_side_coord_atom_positions and surface atom positions
    # Rotate extreme points
    extreme_points = np.array(far_side_coord_atom_positions + [unit_cell.get_positions()[i] for i in surface_atoms])
    extreme_points_centered = extreme_points - center_pos
    rotated_extreme_points = np.dot(extreme_points_centered, R_opt.T) + center_pos
    
    # Compute min and max extremes
    min_extremes = np.min(rotated_extreme_points, axis=0)
    max_extremes = np.max(rotated_extreme_points, axis=0)
    
    # Compute the center of the extreme points
    center_extremes = (min_extremes + max_extremes) / 2
    
    # Compute translation to align center_extremes with cell center
    cell_center = np.array([lattice_constant/2]*3)
    translation = cell_center - center_extremes
    
    # Apply translation to the entire unit cell
    unit_cell.translate(translation)

    # Now adjust the positions to ensure all atoms are within the cell boundaries
    positions = unit_cell.get_positions()

    min_positions = np.min(positions, axis=0)
    max_positions = np.max(positions, axis=0)

    # Initialize additional translation
    additional_translation = np.zeros(3)

    # Adjust along x-axis
    if min_positions[0] < 0:
        additional_translation[0] = -min_positions[0]
    elif max_positions[0] > lattice_constant:
        additional_translation[0] = lattice_constant - max_positions[0]

    # Adjust along y-axis
    if min_positions[1] < 0:
        additional_translation[1] = -min_positions[1]
    elif max_positions[1] > lattice_constant:
        additional_translation[1] = lattice_constant - max_positions[1]

    # Adjust along z-axis
    if min_positions[2] < 0:
        additional_translation[2] = -min_positions[2]
    elif max_positions[2] > lattice_constant:
        additional_translation[2] = lattice_constant - max_positions[2]

    # Apply additional translation
    positions += additional_translation
    unit_cell.set_positions(positions)
    
    metal_formula = metal_center.get_chemical_formula()
    ligand_formula = ligand.get_chemical_formula()
    filename = f"{metal_formula}_{ligand_formula}_cubic_cell.xyz"
    write(filename, unit_cell)
    
    return unit_cell
