# Matbench Discovery data preprocessing

import os
import numpy as np
import torch
from torch_geometric.data import Data
from matbench_discovery.data import load
from pymatgen.core import Structure
import itertools
from tqdm import tqdm
import matterfold as mf


def preprocess_matbench_discovery(path, cutoff: float = 5.0, batch_size: int = 100):
    # Decide which device to use
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    def structures_to_graphs(structures: list, material_ids: list, formulas: list):
        batch_size_local = len(structures)
        max_atoms = max(len(s) for s in structures)

        # Allocate empty tensors directly on the chosen device
        atomic_numbers = torch.zeros((batch_size_local, max_atoms), dtype=torch.long, device=device)
        positions = torch.zeros((batch_size_local, max_atoms, 3), dtype=torch.float32, device=device)
        cells = torch.zeros((batch_size_local, 3, 3), dtype=torch.float32, device=device)
        pbcs = torch.zeros((batch_size_local, 3), dtype=torch.bool, device=device)
        num_atoms = torch.zeros(batch_size_local, dtype=torch.long, device=device)

        # Fill them in with structure data
        for i, structure in enumerate(structures):
            n_atoms = len(structure)
            num_atoms[i] = n_atoms
            atomic_numbers[i, :n_atoms] = torch.tensor(
                [site.specie.number for site in structure],
                dtype=torch.long,
                device=device
            )
            positions[i, :n_atoms] = torch.tensor(
                structure.cart_coords, dtype=torch.float32, device=device
            )
            cells[i] = torch.tensor(structure.lattice.matrix, dtype=torch.float32, device=device)
            pbcs[i] = torch.tensor(structure.lattice.pbc, dtype=torch.bool, device=device)

        # Precompute all PBC offsets on the same device
        pbc_offsets = torch.tensor(
            list(itertools.product([-1, 0, 1], repeat=3)),
            dtype=torch.float32,
            device=device
        )
        # Compute scaled positions
        scaled_positions = torch.matmul(positions, torch.inverse(cells).transpose(1, 2))
        extended_scaled_positions = (
            scaled_positions.unsqueeze(2) + pbc_offsets.unsqueeze(0).unsqueeze(0)
        ).reshape(batch_size_local, -1, 3)
        extended_cart_positions = torch.matmul(extended_scaled_positions, cells)

        diff = extended_cart_positions.unsqueeze(2) - positions.unsqueeze(1)
        dist_matrix = torch.norm(diff, dim=-1)

        graphs = []
        for i in range(batch_size_local):
            # Find neighbors within cutoff distance
            n_atom_27 = num_atoms[i] * 27
            neighbors = torch.nonzero(dist_matrix[i, :n_atom_27, :num_atoms[i]] < cutoff)

            # edge_index: [2, num_edges], separate row and col indices
            edge_index = neighbors[neighbors[:, 0] // 27 != neighbors[:, 1]].T
            # Convert extended atom indices to just 0-based original atom indices
            edge_index[0] = edge_index[0] // 27

            # Build graph data
            x = atomic_numbers[i, :num_atoms[i]].unsqueeze(1).float().to(device)
            pos = positions[i, :num_atoms[i]].to(device)
            cell = cells[i].to(device)
            pbc = pbcs[i].to(device)

            data_graph = Data(
                x=x,
                edge_index=edge_index,
                pos=pos,
                cell=cell,
                pbc=pbc,
                material_id=material_ids[i],
                formula=formulas[i]
            )
            graphs.append(data_graph)

        return graphs

    def process_dataset(structures_data, properties_data, structure_key, energy_key, formula_key, dataset_name):
        print(f"Processing {dataset_name} data...")

        # Intersection of indices to ensure alignment
        common_ids = structures_data.index.intersection(properties_data.index)
        structures_data = structures_data.loc[common_ids]
        properties_data = properties_data.loc[common_ids]
        material_ids = structures_data.index.tolist()

        # Extract structures and formulas
        if dataset_name == 'MP':
            # MP data: 'structure' is nested under `entry`
            structures = [
                Structure.from_dict(row[structure_key]['structure']) 
                for _, row in structures_data.iterrows()
            ]
            formulas = [
                row[structure_key]['composition']
                for _, row in structures_data.iterrows()
            ]
        else:
            # WBM data: 'initial_structure' is directly in the DataFrame row
            structures = [
                Structure.from_dict(row[structure_key]) 
                for _, row in structures_data.iterrows()
            ]
            formulas = properties_data[formula_key].tolist()

        energies = properties_data[energy_key].tolist()

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

                for g, energy in zip(batch_graphs, batch_energies):
                    # Move energy to device as well
                    g.y = torch.tensor([energy], dtype=torch.float, device=device)
                    graphs.append(g)
                    y_values.append(energy)

                pbar.update(1)

        print(f"Processed {dataset_name} structures: {len(graphs)}")
        return graphs, y_values

    # Load MP data
    print("Loading MP data...")
    mp_entries = load("mp_computed_structure_entries", version="1.0.0")
    mp_energies = load("mp_energies", version="1.0.0")
    mp_graphs, mp_y_values = process_dataset(
        mp_entries,
        mp_energies,
        structure_key='entry',
        energy_key='energy_above_hull',
        formula_key='formula_pretty',
        dataset_name='MP'
    )

    # Load WBM data
    print("Loading WBM data...")
    wbm_summary = load("wbm_summary", version="1.0.0")
    wbm_initial_structures = load("wbm_initial_structures", version="1.0.0")
    wbm_graphs, wbm_y_values = process_dataset(
        wbm_initial_structures,
        wbm_summary,
        structure_key='initial_structure',
        energy_key='e_above_hull_mp2020_corrected_ppd_mp',
        formula_key='formula',
        dataset_name='WBM'
    )

    # Save to disk so SOAP function can load/modify in place
    data_save_path = os.path.join(path, 'MBDData.pt')
    print(f"Saving data to {data_save_path}...")
    torch.save({
        'mp_graphs': mp_graphs,
        'mp_y_values': mp_y_values,
        'wbm_graphs': wbm_graphs,
        'wbm_y_values': wbm_y_values
    }, data_save_path)

    # Calculate SOAP descriptors
    print("Calculating SOAP descriptors...")
    mf.soap(
        data_save_path,
        soap_local=True,
        soap_global=True,
        r_cut=5.0,
        n_max=4,
        l_max=4,
        sigma=0.4
    )

    # Reload the data (which now includes soap and global_soap features)
    data = torch.load(data_save_path)
    mp_graphs = data['mp_graphs']
    wbm_graphs = data['wbm_graphs']

    # Combine x=[atomic_number, soap, global_soap] in new graphs
    print("Combining x, soap, and global_soap into new graphs...")

    mp = []
    print("Processing MP graphs with SOAP embeddings...")
    for g in tqdm(mp_graphs, desc="MP Graphs"):
        # Move each field to device and convert to float
        atomic_numbers = g.x.view(-1, 1).float().to(device)
        soap_local = g.soap.float().to(device)
        # global_soap is shape [soap_size], expand to [num_atoms, soap_size]
        soap_global = g.global_soap.unsqueeze(0).repeat(atomic_numbers.size(0), 1).float().to(device)

        new_x = torch.cat([atomic_numbers, soap_local, soap_global], dim=1).to(device)

        new_data = Data(
            x=new_x,
            edge_index=g.edge_index.to(device),
            y=g.y.to(device)
        )
        mp.append(new_data)

    wbm = []
    print("Processing WBM graphs with SOAP embeddings...")
    for g in tqdm(wbm_graphs, desc="WBM Graphs"):
        atomic_numbers = g.x.view(-1, 1).float().to(device)
        soap_local = g.soap.float().to(device)
        soap_global = g.global_soap.unsqueeze(0).repeat(atomic_numbers.size(0), 1).float().to(device)

        new_x = torch.cat([atomic_numbers, soap_local, soap_global], dim=1).to(device)

        new_data = Data(
            x=new_x,
            edge_index=g.edge_index.to(device),
            y=g.y.to(device)
        )
        wbm.append(new_data)

    data['mp'] = mp
    data['wbm'] = wbm

    torch.save(data, data_save_path)
    print(f"Preprocessing complete. Saved to: {data_save_path}.")


if __name__ == "__main__":
    preprocess_matbench_discovery(path="/content/drive/MyDrive/matbench_discovery")

