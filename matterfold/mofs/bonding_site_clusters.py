# ligand bonding-site clustering

import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial import ConvexHull
from ase import Atoms
from ase.geometry import get_distances

def find_surface_bonding_sites(
    ligand: Atoms,
    bonding_sites: dict,
    donor_elements=['N', 'O', 'S', 'P', 'F', 'Cl', 'Se', 'Br', 'I', 'At', 'Ts'],
    min_unbonded_electrons=1.0,
    clustering_distance_threshold=2.5
):
    """
    Identify surface-accessible donor atoms and cluster them into bonding sites.

    Parameters:
    - ligand: ASE Atoms object representing the ligand molecule.
    - bonding_sites: Dictionary containing unbonded electron counts per atom.
    - donor_elements: List of element symbols to consider as potential donors.
    - min_unbonded_electrons: Minimum number of unbonded electrons required.
    - clustering_distance_threshold: Distance threshold for clustering (in Angstroms).

    Returns:
    - clusters: List of clusters, each containing atom indices of a bonding site.
    """
    # Step 1: Filter for potential donor atoms
    donor_atom_indices = []
    donor_positions = []
    for atom_index, info in bonding_sites.items():
        element = info['element']
        non_bonded_electrons = info['non_bonded_electrons']
        if element in donor_elements and non_bonded_electrons >= min_unbonded_electrons:
            donor_atom_indices.append(atom_index - 1)  # Adjusting index to 0-based
            donor_positions.append(ligand.positions[atom_index - 1])

    donor_positions = np.array(donor_positions)

    if len(donor_positions) == 0:
        print("No potential donor atoms found with the specified criteria.")
        return []

    # Step 2: Perform PCA on the ligand positions
    pca = PCA(n_components=3)
    pca.fit(ligand.positions)
    transformed_positions = pca.transform(ligand.positions)
    transformed_donor_positions = pca.transform(donor_positions)

    # Step 3: Determine surface-accessible donor atoms
    # Compute the convex hull of the ligand
    hull = ConvexHull(transformed_positions)
    hull_vertices = hull.vertices

    # Identify donor atoms that are on the convex hull (surface-accessible)
    surface_donor_indices = []
    for idx, pos in zip(donor_atom_indices, transformed_donor_positions):
        # Check if the donor atom index is among the hull vertices
        if idx in hull_vertices:
            surface_donor_indices.append(idx)

    if len(surface_donor_indices) == 0:
        print("No surface-accessible donor atoms found.")
        return []

    surface_donor_positions = ligand.positions[surface_donor_indices]

    # Step 4: Cluster the surface-accessible donor atoms
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=clustering_distance_threshold,
        linkage='single'
    )
    clustering.fit(surface_donor_positions)

    # Organize clusters
    clusters = {}
    for atom_idx, label in zip(surface_donor_indices, clustering.labels_):
        clusters.setdefault(label, []).append(atom_idx)

    # Convert clusters to a list
    bonding_sites_clusters = list(clusters.values())

    return bonding_sites_clusters
