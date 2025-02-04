# TOX21

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
import numpy as np
import torch
import torch.nn.functional as F
from torch.nn import Linear
import sklearn.metrics as metrics
from torch_geometric.datasets import MoleculeNet
from torch_geometric.nn import global_add_pool, GraphConv, GATv2Conv
from torch_geometric.loader import DataLoader


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
path = os.path.join(os.getcwd(), 'data', 'TOX21')
os.makedirs(path, exist_ok=True)

dataset = 'TOX21'
heads = 3
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 500
lr = 0.0005
weight_decay = 0 #1e-5
batch_size = 64

dataset = MoleculeNet(path, name="TOX21").shuffle()

train_dataset = dataset[len(dataset) // 10:]
test_dataset = dataset[:len(dataset) // 10]

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=batch_size)

class CARATE(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels, heads, dropout_gat, dropout_forward):
        super(CARATE, self).__init__()

        self.conv1 = GraphConv(in_channels, hidden_channels)
        self.gat1 = GATv2Conv(hidden_channels, hidden_channels, heads=heads, dropout=dropout_gat, residual=True)
        self.conv2 = GraphConv(hidden_channels * heads, hidden_channels)

        self.fc1 = Linear(hidden_channels, hidden_channels)
        self.fc2 = Linear(hidden_channels, out_channels)

        self.dropout_forward = dropout_forward

    def forward(self, x, edge_index, batch):
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=self.dropout_forward, training=self.training)
        x = F.relu(self.gat1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = global_add_pool(x, batch)
        x = F.relu(self.fc1(x))
        x = F.dropout(x, p=self.dropout_forward, training=self.training)
        x = self.fc2(x)
        return torch.sigmoid(x)

model = CARATE(dataset.num_features, hidden_channels, dataset.num_classes, heads, dropout_gat, dropout_forward).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

def train(epoch):
    model.train()

    if epoch == 51:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.5 * param_group['lr']

    correct = 0
    for data in train_loader:
        data.x = data.x.type(torch.FloatTensor)
        data.y = torch.nan_to_num(data.y.type(torch.FloatTensor))
        data = data.to(device)
        optimizer.zero_grad()
        output_probs = model(data.x, data.edge_index, data.batch)
        output = (output_probs > 0.5).float()
        loss = torch.nn.BCELoss()
        loss = loss(output_probs, data.y)
        loss.backward()
        optimizer.step()
        correct += (output == data.y).float().sum()/dataset.num_classes
    return correct / len(train_loader.dataset)


def test(loader,  epoch, test=False):
    model.eval()

    correct = 0
    auc = np.zeros(dataset.num_classes)
    if test:
      outs =[]
    for data in loader:
        data.x = data.x.type(torch.FloatTensor)
        data.y = torch.nan_to_num(data.y.type(torch.FloatTensor))
        data = data.to(device)
        output_probs = model(data.x, data.edge_index, data.batch)
        output = (output_probs > 0.5).float()
        if test:
          outs.append(output_probs.cpu().detach().numpy())
        correct += (output == data.y).float().sum()/dataset.num_classes
    if test:
      outputs =np.concatenate(outs, axis=0 ).astype(float)
      np.savetxt("Tox21_epoch"+str(epoch)+".csv", outputs)
    return correct / len(loader.dataset)

for epoch in range(1, epochs):
    train_loss = train(epoch)
    train_acc = test(train_loader, epoch)
    test_acc = test(test_loader, test=True, epoch=epoch)
    print('Epoch: {:03d}, Train Loss: {:.7f}, Train Acc: {:.7f}, Test Acc: {:.7f}'.format(epoch, train_loss, train_acc, test_acc))
