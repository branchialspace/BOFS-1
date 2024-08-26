# Matbench-Discovery

import os
import time
import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, GATv2Conv, global_add_pool
from torch_geometric.loader import DataLoader
from matbench_discovery.data import load_df_wbm_with_preds
from matbench_discovery.enums import MbdKey
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from torch_geometric.data import Data
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def structure_to_graph(structure: Structure, cutoff: float = 5.0):
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    
    atomic_numbers = torch.tensor(ase_atoms.get_atomic_numbers(), dtype=torch.long)
    pos = torch.tensor(ase_atoms.get_positions(), dtype=torch.float)
    cell = torch.tensor(ase_atoms.get_cell(), dtype=torch.float)
    pbc = torch.tensor(ase_atoms.get_pbc(), dtype=torch.bool)
    
    edge_index = []
    for i in range(len(atomic_numbers)):
        distances = structure.get_distances(i, range(len(atomic_numbers)), mic=True)
        neighbors = torch.where(torch.tensor(distances) < cutoff)[0]
        edge_index.extend([[i, j] for j in neighbors if i != j])
    
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    
    x = atomic_numbers.unsqueeze(1).float()
    
    return Data(x=x, edge_index=edge_index, pos=pos, cell=cell, pbc=pbc)

df = load_df_wbm_with_preds()

graphs = []
y_values = []
for structure, e_above_hull in zip(df[MbdKey.STRUCTURE], df[MbdKey.E_ABOVE_HULL]):
    graph = structure_to_graph(Structure.from_str(structure, fmt="cif"))
    graph.y = torch.tensor([e_above_hull], dtype=torch.float)
    graphs.append(graph)
    y_values.append(e_above_hull)

train_graphs, test_graphs, train_y, test_y = train_test_split(graphs, y_values, test_size=0.2, random_state=42)

heads = 3
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 50
lr = 0.0005
batch_size = 64

train_loader = DataLoader(train_graphs, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_graphs, batch_size=batch_size)

class CARATE(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, heads, dropout_gat, dropout_forward):
        super().__init__()
        self.conv1 = GraphConv(in_channels, hidden_channels)
        self.gat1 = GATv2Conv(hidden_channels, hidden_channels, heads=heads, dropout=dropout_forward, residual=True)
        self.conv2 = GraphConv(hidden_channels * heads, hidden_channels)
        self.fc1 = torch.nn.Linear(hidden_channels, hidden_channels)
        self.fc2 = torch.nn.Linear(hidden_channels, out_channels)
        self.dropout_gat = dropout_gat
        self.dropout_forward = dropout_forward

    def forward(self, x, edge_index, batch):
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=self.dropout_gat, training=self.training)
        x = F.relu(self.gat1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = global_add_pool(x, batch)
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.dropout_forward, training=self.training)
        x = self.fc2(x)
        return x

model = CARATE(train_graphs[0].x.shape[1], hidden_channels, 1, heads, dropout_gat, dropout_forward).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

def train():
    model.train()
    total_loss = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.batch)
        loss = F.mse_loss(out, data.y)
        loss.backward()
        optimizer.step()
        total_loss += float(loss) * data.num_graphs
    return total_loss / len(train_loader.dataset)

@torch.no_grad()
def test(loader):
    model.eval()
    total_mae = 0
    all_preds = []
    all_true = []
    for data in loader:
        data = data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        total_mae += F.l1_loss(out, data.y).item() * data.num_graphs
        all_preds.extend(out.cpu().numpy())
        all_true.extend(data.y.cpu().numpy())
    
    mae = total_mae / len(loader.dataset)
    
    all_preds = np.array(all_preds).flatten()
    all_true = np.array(all_true).flatten()
    pred_stable = (all_preds < 0.1).astype(int)
    true_stable = (all_true < 0.1).astype(int)
    f1 = f1_score(true_stable, pred_stable)
    
    return mae, f1

best_val_error = float('inf')
test_mae = test_f1 = 0

for epoch in range(1, epochs + 1):
    start = time.time()
    loss = train()
    train_mae, train_f1 = test(train_loader)
    val_mae, val_f1 = test(test_loader)
    
    if val_mae < best_val_error:
        best_val_error = val_mae
        test_mae = val_mae
        test_f1 = val_f1

    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train MAE: {train_mae:.4f}, '
          f'Train F1: {train_f1:.4f}, Val MAE: {val_mae:.4f}, Val F1: {val_f1:.4f}')

    if epoch % 10 == 0:
        torch.save(model.state_dict(), f'model_epoch_{epoch}.pt')

print(f'Final Test MAE: {test_mae:.4f}, Test F1: {test_f1:.4f}')
