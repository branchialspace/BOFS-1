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

    def structures_to_graphs(structures: list, material_ids: list, formulas: list):
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
            graphs.append(Data(x=x, edge_index=edge_index, pos=pos, cell=cell, pbc=pbc, 
                               material_id=material_ids[i], formula=formulas[i]))

        return graphs

    def process_dataset(structures_data, properties_data, structure_key, energy_key, formula_key, dataset_name):
        print(f"Processing {dataset_name} data...")
        
        # Ensure data alignment
        common_ids = structures_data.index.intersection(properties_data.index)
        structures_data = structures_data.loc[common_ids]
        properties_data = properties_data.loc[common_ids]

        structures = [Structure.from_dict(row[structure_key]) for _, row in structures_data.iterrows()]
        energies = properties_data[energy_key].tolist()
        material_ids = properties_data.index.tolist()
        formulas = properties_data[formula_key].tolist()

        print(f"Total {dataset_name} structures: {len(structures)}")

        graphs = []
        y_values = []

        total_batches = (len(structures) + batch_size - 1) // batch_size
        with tqdm(total=total_batches, desc=f"{dataset_name} Batches") as pbar:
            for i in range(0, len(structures), batch_size):
                batch_structures = structures[i:i+batch_size]
                batch_material_ids = material_ids[i:i+batch_size]
                batch_formulas = formulas[i:i+batch_size]
                batch_graphs = structures_to_graphs(batch_structures, batch_material_ids, batch_formulas)
                batch_energies = energies[i:i+batch_size]

                for graph, energy in zip(batch_graphs, batch_energies):
                    graph.y = torch.tensor([energy], dtype=torch.float)
                    graphs.append(graph)
                    y_values.append(energy)

                pbar.update(1)

        print(f"Processed {dataset_name} structures: {len(graphs)}")
        return graphs, y_values

    # Process MP data
    mp_entries = load("mp_computed_structure_entries", version="1.0.0")
    mp_energies = load("mp_energies", version="1.0.0")
    mp_graphs, mp_y_values = process_dataset(
        mp_entries, mp_energies, 
        structure_key='entry', 
        energy_key='energy_above_hull', 
        formula_key='formula_pretty',
        dataset_name='MP'
    )

    # Process WBM data
    wbm_summary = load("wbm_summary", version="1.0.0")
    wbm_initial_structures = load("wbm_initial_structures", version="1.0.0")
    
    # Prepare WBM data to match MP data structure
    wbm_structures = wbm_initial_structures.set_index('formula_from_cse')
    wbm_structures = wbm_structures.rename(columns={'initial_structure': 'structure'})
    wbm_summary = wbm_summary.set_index('formula')
    
    wbm_graphs, wbm_y_values = process_dataset(
        wbm_structures, wbm_summary, 
        structure_key='structure', 
        energy_key='e_above_hull_mp2020_corrected_ppd_mp', 
        formula_key='formula',
        dataset_name='WBM'
    )

    print(f"Saving processed data to {path}...")
    torch.save({
        'mp_graphs': mp_graphs,
        'mp_y_values': mp_y_values,
        'wbm_graphs': wbm_graphs,
        'wbm_y_values': wbm_y_values
    }, os.path.join(path, 'MBDData.pt'))

    print(f"Preprocessing complete. Data saved as 'MBDData.pt' (cutoff distance: {cutoff})")

if __name__ == "__main__":
    preprocess_matbench_discovery_data(path='/content/drive/MyDrive/matbench_discovery/')
