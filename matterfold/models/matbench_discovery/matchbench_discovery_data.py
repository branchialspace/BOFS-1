# Matbench Discovery data preprocessing

import os
import numpy as np
import torch
from torch_geometric.data import Data
from matbench_discovery.data import load
from pymatgen.core import Structure
import itertools
from tqdm import tqdm

def preprocess_matbench_discovery_data(path, cutoff: float = 5.0, batch_size: int = 100):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    def structures_to_graphs(structures: list):
        batch_size = len(structures)
        max_atoms = max(len(s) for s in structures)

        atomic_numbers = torch.zeros((batch_size, max_atoms), dtype=torch.long, device=device)
        positions = torch.zeros((batch_size, max_atoms, 3), dtype=torch.float32, device=device)
        cells = torch.zeros((batch_size, 3, 3), dtype=torch.float32, device=device)
        pbcs = torch.zeros((batch_size, 3), dtype=torch.bool, device=device)
        num_atoms = torch.zeros(batch_size, dtype=torch.long, device=device)

        for i, structure in enumerate(structures):
            n_atoms = len(structure)
            num_atoms[i] = n_atoms
            atomic_numbers[i, :n_atoms] = torch.tensor([site.specie.number for site in structure], dtype=torch.long)
            positions[i, :n_atoms] = torch.tensor(structure.cart_coords, dtype=torch.float32)
            cells[i] = torch.tensor(structure.lattice.matrix, dtype=torch.float32)
            pbcs[i] = torch.tensor(structure.lattice.pbc, dtype=torch.bool)

        pbc_offsets = torch.tensor(list(itertools.product([-1, 0, 1], repeat=3)), dtype=torch.float32, device=device)
        scaled_positions = torch.matmul(positions, torch.inverse(cells).transpose(1, 2))
        extended_scaled_positions = (scaled_positions.unsqueeze(2) + pbc_offsets.unsqueeze(0).unsqueeze(0)).reshape(batch_size, -1, 3)
        extended_cart_positions = torch.matmul(extended_scaled_positions, cells)

        diff = extended_cart_positions.unsqueeze(2) - positions.unsqueeze(1)
        dist_matrix = torch.norm(diff, dim=-1)

        graphs = []
        for i in range(batch_size):
            neighbors = torch.nonzero(dist_matrix[i, :num_atoms[i]*27, :num_atoms[i]] < cutoff)
            edge_index = neighbors[neighbors[:, 0] // 27 != neighbors[:, 1]].T
            edge_index[0] = edge_index[0] // 27
            x = atomic_numbers[i, :num_atoms[i]].unsqueeze(1).float()
            pos = positions[i, :num_atoms[i]]
            cell = cells[i]
            pbc = pbcs[i]
            graphs.append(Data(x=x, edge_index=edge_index, pos=pos, cell=cell, pbc=pbc))

        return graphs

    mp_entries = load("mp_computed_structure_entries", version="1.0.0")
    mp_energies = load("mp_energies", version="1.0.0")

    mp_graphs = []
    mp_y_values = []
    print("Processing MP data...")
    structures = [Structure.from_dict(row['entry']['structure']) for _, row in mp_entries.iterrows()]
    energies = mp_energies['energy_above_hull'].tolist()

    total_batches = (len(structures) + batch_size - 1) // batch_size
    with tqdm(total=total_batches, desc="MP Batches") as pbar:
        for i in range(0, len(structures), batch_size):
            batch_structures = structures[i:i+batch_size]
            batch_graphs = structures_to_graphs(batch_structures)
            batch_energies = energies[i:i+batch_size]

            for graph, energy in zip(batch_graphs, batch_energies):
                graph.y = torch.tensor([energy], dtype=torch.float)
                mp_graphs.append(graph)
                mp_y_values.append(energy)

            pbar.update(1)

    wbm_summary = load("wbm_summary", version="1.0.0")
    wbm_initial_structures = load("wbm_initial_structures", version="1.0.0")

    wbm_graphs = []
    wbm_y_values = []
    print("Processing WBM data...")
    structures = []
    energies = []

    wbm_initial_structures_dict = wbm_initial_structures.set_index('formula')['initial_structure'].to_dict()

    def process_row(row):
        initial_structure = wbm_initial_structures_dict.get(row['formula'])
        if initial_structure is not None:
            return Structure.from_dict(initial_structure), row['e_above_hull_mp2020_corrected_ppd_mp']
        return None, None

    tqdm.pandas(desc="Processing WBM rows")
    results = wbm_summary.progress_apply(process_row, axis=1)
    structures, energies = zip(*results)

    structures = [s for s in structures if s is not None]
    energies = [e for e in energies if e is not None]

    total_batches = (len(structures) + batch_size - 1) // batch_size
    with tqdm(total=total_batches, desc="WBM Batches") as pbar:
        for i in range(0, len(structures), batch_size):
            batch_structures = structures[i:i+batch_size]
            batch_graphs = structures_to_graphs(batch_structures)
            batch_energies = energies[i:i+batch_size]

            for graph, energy in zip(batch_graphs, batch_energies):
                graph.y = torch.tensor([energy], dtype=torch.float)
                wbm_graphs.append(graph)
                wbm_y_values.append(energy)

            pbar.update(1)

    torch.save({
        'mp_graphs': mp_graphs,
        'mp_y_values': mp_y_values,
        'wbm_graphs': wbm_graphs,
        'wbm_y_values': wbm_y_values
    }, os.path.join(path, 'MBDData.pt'))

    print(f"Preprocessing complete. Data saved as 'MBDData.pt' (cutoff distance: {cutoff})")

if __name__ == "__main__":
    preprocess_matbench_discovery_data(path='/content/drive/MyDrive/matbench_discovery/')
