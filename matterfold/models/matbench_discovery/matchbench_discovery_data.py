# Matbench Discovery data preprocessing

import os
import numpy as np
import torch
import torch.multiprocessing as mp
from torch.utils.data import Dataset, DataLoader
from torch_geometric.data import Data, Batch
from matbench_discovery.data import load
from pymatgen.core import Structure
import itertools


class MatbenchDataset(Dataset):
    def __init__(self, entries, energies, is_mp=True):
        self.entries = entries
        self.energies = energies
        self.is_mp = is_mp
        self.device = torch.device('cpu')

        if self.is_mp:
            self.entries = self.entries.reset_index(drop=True)
            self.energies = self.energies.set_index(self.entries.index)

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, idx):
        if self.is_mp:
            entry = self.entries.iloc[idx]['entry']
            structure = Structure.from_dict(entry['structure'])
            energy = self.energies.loc[self.entries.index[idx], 'energy_above_hull']
        else:
            row = self.entries.iloc[idx]
            cse_entry = self.energies[self.energies['formula_from_cse'] == row['formula']].iloc[0]
            structure = Structure.from_dict(cse_entry['computed_structure_entry']['structure'])
            energy = row['e_above_hull_mp2020_corrected_ppd_mp']

        return structure, torch.tensor(energy, dtype=torch.float32)

def structure_to_graph(structure: Structure, cutoff: float = 5.0):
    device = torch.device('cpu')

    atomic_numbers = torch.tensor([site.specie.number for site in structure], dtype=torch.long, device=device)
    pos = torch.tensor(structure.cart_coords, dtype=torch.float32, device=device)
    cell = torch.tensor(structure.lattice.matrix, dtype=torch.float32, device=device)
    pbc = torch.tensor(structure.lattice.pbc, dtype=torch.bool, device=device)

    lattice = structure.lattice
    scaled_positions = torch.tensor(lattice.get_fractional_coords(structure.cart_coords), dtype=torch.float32, device=device)
    pbc_offsets = torch.tensor(list(itertools.product([-1, 0, 1], repeat=3)), dtype=torch.float32, device=device)

    extended_scaled_positions = (scaled_positions.unsqueeze(1) + pbc_offsets.unsqueeze(0)).reshape(-1, 3)
    extended_cart_positions = torch.matmul(extended_scaled_positions, cell.T)

    diff = extended_cart_positions.unsqueeze(1) - pos.unsqueeze(0)
    dist_matrix = torch.norm(diff, dim=-1)
    neighbors = torch.nonzero(dist_matrix < cutoff)
    edge_index = neighbors[neighbors[:, 0] // len(pbc_offsets) != neighbors[:, 1]].T
    edge_index[0] = edge_index[0] // len(pbc_offsets)

    x = atomic_numbers.unsqueeze(1).float()

    return Data(x=x, edge_index=edge_index, pos=pos, cell=cell, pbc=pbc)

def collate_fn(batch):
    structures, energies = zip(*batch)
    graphs = [structure_to_graph(structure) for structure in structures]
    batch_data = Batch.from_data_list(graphs)
    batch_data.y = torch.stack(energies).view(-1, 1)
    return batch_data

def process_dataset(dataset, loader, device):
    graphs = []
    y_values = []

    for i, batch in enumerate(loader):
        batch = batch.to(device)
        graphs.extend(batch.to_data_list())
        y_values.extend(batch.y.cpu().numpy().flatten())
        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1} batches out of {len(loader)}")

    return graphs, y_values

def preprocess_matbench_discovery_data(cutoff: float = 5.0, batch_size: int = 32, num_workers: int = 4):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    mp_entries = load("mp_computed_structure_entries", version="1.0.0")
    mp_energies = load("mp_energies", version="1.0.0")

    mp_dataset = MatbenchDataset(mp_entries, mp_energies, is_mp=True)
    mp_loader = DataLoader(mp_dataset, batch_size=batch_size, num_workers=num_workers, collate_fn=collate_fn)

    print("Processing MP data...")
    mp_graphs, mp_y_values = process_dataset(mp_dataset, mp_loader, device)

    wbm_summary = load("wbm_summary", version="1.0.0")
    wbm_cses = load("wbm_computed_structure_entries", version="1.0.0")

    wbm_dataset = MatbenchDataset(wbm_summary, wbm_cses, is_mp=False)
    wbm_loader = DataLoader(wbm_dataset, batch_size=batch_size, num_workers=num_workers, collate_fn=collate_fn)

    print("Processing WBM data...")
    wbm_graphs, wbm_y_values = process_dataset(wbm_dataset, wbm_loader, device)

    torch.save({
        'mp_graphs': mp_graphs,
        'mp_y_values': mp_y_values,
        'wbm_graphs': wbm_graphs,
        'wbm_y_values': wbm_y_values
    }, 'preprocessed_data.pt')

    print(f"Preprocessing complete. Data saved to 'preprocessed_data.pt' (cutoff distance: {cutoff})")

if __name__ == "__main__":
    mp.set_start_method('spawn')
    preprocess_matbench_discovery_data()
