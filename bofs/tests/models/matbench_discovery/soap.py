# Matbench Discovery SOAP positional encodings

import torch
import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP
from tqdm import tqdm


def soap(data_path, soap_local=True, soap_global=True, r_cut=5.0, n_max=4, l_max=4, sigma=0.4):
    data = torch.load(data_path)

    # Initialize SOAP descriptors
    soap_descriptor = None
    global_soap_descriptor = None

    if soap_local:
        soap_descriptor = SOAP(species=["H"], periodic=True, r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma)

    if soap_global:
        global_soap_descriptor = SOAP(species=["H"], periodic=True, r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma, average="inner")

    # Process MP graphs
    for graph in tqdm(data['mp_graphs'], desc="Processing MP graphs"):
        process_graph(graph, soap=soap_descriptor, average_soap=global_soap_descriptor)

    # Process WBM graphs
    for graph in tqdm(data['wbm_graphs'], desc="Processing WBM graphs"):
        process_graph(graph, soap=soap_descriptor, average_soap=global_soap_descriptor)

    # Save data once
    torch.save(data, data_path)
    print(f"SOAP descriptors calculated and saved to {data_path}.")

def process_graph(graph, soap=None, average_soap=None):
    positions = graph.pos.cpu().numpy()
    cell = graph.cell.cpu().numpy() if isinstance(graph.cell, torch.Tensor) else graph.cell
    pbc = graph.pbc.cpu().numpy() if isinstance(graph.pbc, torch.Tensor) else graph.pbc

    # Ensure pbc is boolean
    if isinstance(pbc, np.ndarray) and pbc.dtype != bool:
        pbc = pbc.astype(bool)

    # Create ASE system
    system = Atoms(
        numbers=np.ones(len(positions), dtype=int),
        positions=positions,
        cell=cell,
        pbc=pbc
    )

    # Compute local SOAP
    if soap is not None:
        soap_descriptors = soap.create(system)
        graph.soap = torch.tensor(soap_descriptors, dtype=torch.float32)

    # Compute global SOAP
    if average_soap is not None:
        global_soap_descriptor = average_soap.create(system)
        graph.global_soap = torch.tensor(global_soap_descriptor, dtype=torch.float32)

if __name__ == "__main__":

    data_path = '/content/drive/MyDrive/matbench_discovery/MBDData.pt'
    soap(data_path, soap_local=True, soap_global=True, r_cut=5.0, n_max=4, l_max=4, sigma=0.4)
