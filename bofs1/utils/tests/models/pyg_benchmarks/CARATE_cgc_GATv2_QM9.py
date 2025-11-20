#QM9

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
import torch_geometric
from torch_geometric.datasets import QM9
from torch_geometric.logging import init_wandb, log
from torch_geometric.nn import GraphConv, GATv2Conv, global_add_pool
from torch_geometric.loader import DataLoader


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
path = os.path.join(os.getcwd(), 'data', 'QM9')
os.makedirs(path, exist_ok=True)

dataset = 'QM9'
heads = 1
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 5000
lr = 0.000005
batch_size = 32

init_wandb(name=f'CARATE-{dataset}', heads=heads, epochs=epochs, hidden_channels=hidden_channels, lr=lr, device=device)

full_dataset = QM9(path).shuffle()
train_dataset = full_dataset[len(full_dataset) // 10:]
test_dataset = full_dataset[:len(full_dataset) // 10]

y = full_dataset.data.y.numpy()
factor = np.zeros((1, y.shape[1]))
for i in range(y.shape[1]):
    norm = np.linalg.norm(y[:, i], ord=2)
    factor[0, i] = norm

print(f"Normalization factor: {factor}")

class Normalize(object):
    def __call__(self, data):
        data.x = data.x.float()
        data.y = (data.y / torch.from_numpy(factor)).float()
        return data

transform = Normalize()

train_dataset.transform = transform
test_dataset.transform = transform

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=batch_size)

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

model = CARATE(full_dataset.num_features, hidden_channels, full_dataset.num_classes, heads, dropout_gat, dropout_forward).to(device)
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
    mae = 0
    mse = 0
    for data in loader:
        data = data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        mae += F.l1_loss(out, data.y).item() * data.num_graphs
        mse += F.mse_loss(out, data.y).item() * data.num_graphs
    return mae / len(loader.dataset), mse / len(loader.dataset)

times = []
best_val_error = float('inf')
test_mae = float('inf')
test_mse = float('inf')

for epoch in range(1, epochs + 1):
    start = time.time()
    loss = train()
    train_mae, train_mse = test(train_loader)
    test_mae, test_mse = test(test_loader)
    log(Epoch=epoch, Loss=loss, Train_MAE=train_mae, Train_MSE=train_mse, Test_MAE=test_mae, Test_MSE=test_mse)
    times.append(time.time() - start)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.7f}, Train MAE: {train_mae:.7f}, Test MAE: {test_mae:.7f}')

    if epoch % 51 == 0:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.5 * param_group['lr']

print(f'Median time per epoch: {torch.tensor(times).median():.4f}s')
print(f'Final Test MAE: {test_mae * factor.mean():.4f}, MSE: {test_mse * (factor.mean()**2):.4f}')
