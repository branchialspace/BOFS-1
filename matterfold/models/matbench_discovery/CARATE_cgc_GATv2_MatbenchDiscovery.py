# Matbench Discovery

"""
Implementation of CARATE (Chemistry Attention Representation through Encoders)
Based on the following works:

Academic Paper:
Kleber, J. M. (2024). Introducing CARATE: Finally Speaking Chemistry Through Learning 
Hidden Wave-Function Representations on Graph-Based Attention and Convolutional Neural Networks. 
Research & Reviews: Journal of Chemistry, 13, 1-10.
DOI: 10.4172/2319-9849.13.3.001

Original Implementation:
Kleber, J. M. (2021-2024). CARATE - Chemistry Attention Representation through Encoders.
Repository: https://codeberg.org/sail.black/carate
"""

import os
import time
import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.nn import GraphConv, GATv2Conv, global_add_pool
from torch_geometric.loader import DataLoader
from sklearn.metrics import f1_score

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

heads = 3
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 50
lr = 0.0005
batch_size = 64

# Processed Matbench Discovery training data with energy above hull graph-level regression targets
# !gdown 1--ibJyQ3SNnTSHjdcSocy-QxOOO7Dsa5
# !unzip MBDData.zip -d ./MBDData/ > /dev/null 2>&1

# Load preprocessed data
preprocessed_data = torch.load('/content/MBDData/MBDData.pt')
mp_graphs = preprocessed_data['mp_graphs']
mp_y_values = preprocessed_data['mp_y_values']
wbm_graphs = preprocessed_data['wbm_graphs']
wbm_y_values = preprocessed_data['wbm_y_values']

# Use MP data for training
train_loader = DataLoader(mp_graphs, batch_size=batch_size, shuffle=True)

# Use WBM data for testing
test_loader = DataLoader(wbm_graphs, batch_size=batch_size)

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

model = CARATE(mp_graphs[0].x.shape[1], hidden_channels, 1, heads, dropout_gat, dropout_forward).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

def train():
    model.train()
    total_loss = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.batch)
        loss = F.mse_loss(out.squeeze(), data.y)
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
        total_mae += F.l1_loss(out.squeeze(), data.y).item() * data.num_graphs
        all_preds.extend(out.cpu().numpy())
        all_true.extend(data.y.cpu().numpy())

    mae = total_mae / len(loader.dataset)

    all_preds = np.array(all_preds).flatten()
    all_true = np.array(all_true).flatten()
    pred_stable = (all_preds <= 0).astype(int)
    true_stable = (all_true <= 0).astype(int)
    f1 = f1_score(true_stable, pred_stable)

    return mae, f1

best_val_error = float('inf')
test_mae = test_f1 = 0

for epoch in range(1, epochs + 1):
    start = time.time()
    loss = train()
    train_mae, train_f1 = test(train_loader)
    test_mae, test_f1 = test(test_loader)

    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train MAE: {train_mae:.4f}, '
          f'Train F1: {train_f1:.4f}, Test MAE: {test_mae:.4f}, Test F1: {test_f1:.4f}')

    if epoch % 10 == 0:
        torch.save(model.state_dict(), f'model_epoch_{epoch}.pt')

print(f'Final Test MAE: {test_mae:.4f}, Test F1: {test_f1:.4f}')

# Calculate Discovery Acceleration Factor (DAF)
def calculate_daf(predictions, true_values, k=10000):
    sorted_indices = np.argsort(predictions)
    top_k_true = np.sum(true_values[sorted_indices[:k]] <= 0)
    random_k_true = np.sum(true_values <= 0) * (k / len(true_values))
    return top_k_true / random_k_true

model.eval()
all_preds = []
all_true = []
for data in test_loader:
    data = data.to(device)
    out = model(data.x, data.edge_index, data.batch)
    all_preds.extend(out.cpu().numpy())
    all_true.extend(data.y.cpu().numpy())

all_preds = np.array(all_preds).flatten()
all_true = np.array(all_true).flatten()

daf = calculate_daf(all_preds, all_true)
print(f'Discovery Acceleration Factor (DAF): {daf:.2f}')
