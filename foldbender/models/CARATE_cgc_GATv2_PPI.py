# PPI

import os
import numpy as np
import torch
import torch.nn.functional as F
from torch.nn import Linear
import sklearn.metrics as metrics
from torch_geometric.datasets import PPI
from torch_geometric.nn import global_add_pool, GATv2Conv, GraphConv
from torch_geometric.loader import DataLoader


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
path = os.path.join(os.getcwd(), 'data', 'PPI')
os.makedirs(path, exist_ok=True)

heads = 3
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 5000
lr = 0.0005
weight_decay = 0  # 1e-5
batch_size = 1

train_dataset = PPI(path, split='train')
val_dataset = PPI(path, split='val')
test_dataset = PPI(path, split='test')

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size)
test_loader = DataLoader(test_dataset, batch_size=batch_size)

train_num_nodes = sum(data.num_nodes for data in train_dataset)

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
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.dropout_forward, training=self.training)
        x = self.fc2(x)
        return x

model = CARATE(train_dataset.num_features, hidden_channels, train_dataset.num_classes, heads, dropout_gat, dropout_forward).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

def train(epoch):
    model.train()
    
    if epoch == 51:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.5 * param_group['lr']
    
    total_loss = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.batch)
        loss = F.binary_cross_entropy_with_logits(out, data.y)
        total_loss += loss.item() * data.num_nodes
        loss.backward()
        optimizer.step()
    return total_loss / train_num_nodes

@torch.no_grad()
def test(loader):
    model.eval()
    
    ys, preds = [], []
    for data in loader:
        data = data.to(device)
        out = model(data.x, data.edge_index, data.batch)
        ys.append(data.y)
        preds.append((out > 0).float())
    
    y, pred = torch.cat(ys, dim=0).cpu().numpy(), torch.cat(preds, dim=0).cpu().numpy()
    return metrics.f1_score(y, pred, average='micro') if pred.sum() > 0 else 0

for epoch in range(1, epochs + 1):
    loss = train(epoch)
    val_f1 = test(val_loader)
    test_f1 = test(test_loader)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Val F1: {val_f1:.4f}, Test F1: {test_f1:.4f}')
