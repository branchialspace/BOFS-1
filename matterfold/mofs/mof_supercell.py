# MOF cubic lattice supercell

from ase import Atoms
from ase.io import write
from ase.build.supercells import find_optimal_cell_shape, make_supercell
import numpy as np


def mof_supercell(unit_cell: Atoms, target_size: int = 8, target_shape: str = 'sc', wrap: bool = True) -> Atoms:
    """
    Creates an extended lattice structure from a MOF unit cell.
    
    Parameters:
    - unit_cell: ASE Atoms object of the MOF primitive cubic unit cell
    - target_size: Target number of unit cells in the supercell (default: 8 for 2x2x2)
    - target_shape: Desired supercell shape, 'sc' for simple cubic or 'fcc' for face-centered cubic (default: 'sc')
    - wrap: Whether to wrap atoms at the boundaries (default: True)
    
    Returns:
    - extended_lattice: ASE Atoms object of the extended MOF lattice
    """
    # Find the optimal transformation matrix for the desired supercell
    P = find_optimal_cell_shape(
        cell=unit_cell.cell,
        target_size=target_size,
        target_shape=target_shape,
        verbose=True
    )
    
    if P is None:
        raise ValueError(f"Could not find suitable transformation matrix for target_size={target_size}")
    
    # Create the supercell using the transformation matrix
    extended_lattice = make_supercell(
        prim=unit_cell,
        P=P,
        wrap=wrap,
        order="cell-major"  # Keep atoms from same unit cell together
    )
    
    formula = unit_cell.get_chemical_formula()
    filename = f"{formula}_extended_{target_size}_{target_shape}.xyz"
    write(filename, extended_lattice)
    
    return extended_lattice
