# Matbench Discovery SOAP positional encodings

import torch
import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP
from tqdm import tqdm


def soap_local(data_path, r_cut=5.0, n_max=4, l_max=4, sigma=0.4):
    data = torch.load(data_path)
    
    soap = SOAP(species=["H"], periodic=True, r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma)
    for graph in tqdm(data['mp_graphs'], desc="Processing MP graphs"):
        process_graph(graph, soap)
    for graph in tqdm(data['wbm_graphs'], desc="Processing WBM graphs"):
        process_graph(graph, soap)
    torch.save(data, data_path)
    print(f"SOAP descriptors calculated and saved to {data_path}.")

def soap_global(data_path, r_cut=5.0, n_max=4, l_max=4, sigma=0.4):
    data = torch.load(data_path)
    
    average_soap = SOAP(species=["H"], periodic=True, r_cut=r_cut, n_max=n_max, l_max=l_max, sigma=sigma, average="inner")
    for graph in tqdm(data['mp_graphs'], desc="Processing MP graphs (global)"):
        process_graph_global(graph, average_soap)
    for graph in tqdm(data['wbm_graphs'], desc="Processing WBM graphs (global)"):
        process_graph_global(graph, average_soap)
    torch.save(data, data_path)
    print(f"Global SOAP descriptors calculated and saved to {data_path}.")

def process_graph(graph, soap):
    positions = graph.pos.numpy()
    num_atoms = positions.shape[0]
    system = Atoms(numbers=np.ones(num_atoms), positions=positions, cell=graph.cell, pbc=graph.pbc)
    soap_descriptors = soap.create(system)
    
    graph.soap = torch.tensor(soap_descriptors, dtype=torch.float32)

def process_graph_global(graph, average_soap):
    positions = graph.pos.numpy()
    num_atoms = positions.shape[0]
    system = Atoms(numbers=np.ones(num_atoms), positions=positions, cell=graph.cell, pbc=graph.pbc)
    global_soap_descriptor = average_soap.create(system)
    
    graph.global_soap = torch.tensor(global_soap_descriptor, dtype=torch.float32)

if __name__ == "__main__":
    data_path = '/content/drive/MyDrive/matbench_discovery/MBDData.pt'
    soap_local(data_path, r_cut=5.0, n_max=4, l_max=4, sigma=0.4)
    soap_global(data_path, r_cut=5.0, n_max=4, l_max=4, sigma=0.4)
