# Matbench Discovery data preprocessing

import os
import numpy as np
import torch
from torch_geometric.data import Data
from matbench_discovery.data import load
from pymatgen.core import Structure
import itertools
from tqdm import tqdm


def preprocess_matbench_discovery_data(cutoff: float = 5.0):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    def structure_to_graph(structure: Structure):
        # Get atomic numbers
        atomic_numbers = torch.tensor([site.specie.number for site in structure], dtype=torch.long)

        # Get positions, cell and pbc
        pos = torch.tensor(structure.cart_coords, dtype=torch.float32).to(device)
        cell = torch.tensor(structure.lattice.matrix, dtype=torch.float32).to(device)
        pbc = torch.tensor(structure.lattice.pbc, dtype=torch.bool).to(device)

        # Generate periodic images
        lattice = structure.lattice
        scaled_positions = torch.tensor(lattice.get_fractional_coords(structure.cart_coords), dtype=torch.float32).to(device)
        pbc_offsets = torch.tensor(list(itertools.product([-1, 0, 1], repeat=3)), dtype=torch.float32).to(device)

        extended_scaled_positions = (scaled_positions.unsqueeze(1) + pbc_offsets.unsqueeze(0)).reshape(-1, 3)
        extended_cart_positions = torch.matmul(extended_scaled_positions, cell.T)

        # Compute pairwise distances
        diff = extended_cart_positions.unsqueeze(1) - pos.unsqueeze(0)
        dist_matrix = torch.norm(diff, dim=-1)
        neighbors = torch.nonzero(dist_matrix < cutoff)
        edge_index = neighbors[neighbors[:, 0] // len(pbc_offsets) != neighbors[:, 1]].T
        edge_index[0] = edge_index[0] // len(pbc_offsets)

        x = atomic_numbers.unsqueeze(1).float().to(device)

        return Data(x=x, edge_index=edge_index, pos=pos, cell=cell, pbc=pbc)

    # Load the MP training data
    mp_entries = load("mp_computed_structure_entries", version="1.0.0")
    mp_energies = load("mp_energies", version="1.0.0")

    # Process MP data
    mp_graphs = []
    mp_y_values = []

    print("Processing MP data...")
    for idx, row in tqdm(mp_entries.iterrows(), total=len(mp_entries)):
        entry = row['entry']
        structure = Structure.from_dict(entry['structure'])
        graph = structure_to_graph(structure)

        # Get the corresponding energy
        energy = mp_energies.loc[idx, 'energy_above_hull']

        graph.y = torch.tensor([energy], dtype=torch.float)
        mp_graphs.append(graph)
        mp_y_values.append(energy)

    # Load the WBM test data
    wbm_summary = load("wbm_summary", version="1.0.0")
    wbm_cses = load("wbm_computed_structure_entries", version="1.0.0")

    # Process WBM data
    wbm_graphs = []
    wbm_y_values = []

    print("Processing WBM data...")
    for idx, row in tqdm(wbm_summary.iterrows(), total=len(wbm_summary)):
        cse_entry = wbm_cses[wbm_cses['formula_from_cse'] == row['formula']].iloc[0]
        structure = Structure.from_dict(cse_entry['computed_structure_entry']['structure'])

        graph = structure_to_graph(structure)
        graph.y = torch.tensor([row['e_above_hull_mp2020_corrected_ppd_mp']], dtype=torch.float)
        wbm_graphs.append(graph)
        wbm_y_values.append(row['e_above_hull_mp2020_corrected_ppd_mp'])

    # Save preprocessed data
    torch.save({
        'mp_graphs': mp_graphs,
        'mp_y_values': mp_y_values,
        'wbm_graphs': wbm_graphs,
        'wbm_y_values': wbm_y_values
    }, 'preprocessed_data.pt')

    print(f"Preprocessing complete. Data saved to 'preprocessed_data.pt' (cutoff distance: {cutoff})")

