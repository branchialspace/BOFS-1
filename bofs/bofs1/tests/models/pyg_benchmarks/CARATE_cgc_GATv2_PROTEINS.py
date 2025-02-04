# PROTEINS

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
from torch_geometric.datasets import TUDataset
from torch_geometric.nn import global_add_pool, GraphConv, GATv2Conv
from torch_geometric.loader import DataLoader


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
path = os.path.join(os.getcwd(), 'data', 'PROTEINS')
os.makedirs(path, exist_ok=True)

dataset = 'PROTEINS_full'
heads = 3
hidden_channels = 364
dropout_gat = 0.6
dropout_forward = 0.5
epochs = 5000
lr = 0.0005
batch_size = 64

dataset = TUDataset(path, name="PROTEINS_full").shuffle()

train_dataset = dataset[len(dataset) // 20:]
test_dataset = dataset[:len(dataset) // 20]

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
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

def train(epoch):
    model.train()

    if epoch == 51:
        for param_group in optimizer.param_groups:
            param_group['lr'] = 0.5 * param_group['lr']

    correct = 0
    for data in train_loader:
        data.x = data.x.type(torch.FloatTensor)
        data.y = F.one_hot(data.y, num_classes=dataset.num_classes).type(torch.FloatTensor)
        data = data.to(device)
        optimizer.zero_grad()
        output_probs = model(data.x, data.edge_index, data.batch)
        output = (output_probs > 0.5).float()
        loss = torch.nn.BCELoss()
        loss = loss(output_probs, data.y)
        loss.backward()
        optimizer.step()
        correct += (output == data.y).float().sum() / dataset.num_classes
    return correct / len(train_loader.dataset), loss

def test(loader, epoch, test=False):
    model.eval()

    correct = 0
    if test:
        outs = []
    for data in loader:
        data.x = data.x.type(torch.FloatTensor)
        data = data.to(device)
        output_probs = model(data.x, data.edge_index, data.batch)
        output = (output_probs > 0.5).float()
        y_one_hot = F.one_hot(data.y, num_classes=dataset.num_classes).float().to(device)
        correct += (output == y_one_hot).float().sum() / dataset.num_classes
        if test:
            outs.append(output.cpu().detach().numpy())
    if test:
        outputs = np.concatenate(outs, axis=0).astype(float)
        np.savetxt(f"PROTEINS_epoch{epoch}.csv", outputs)
    return correct / len(loader.dataset)

for epoch in range(1, epochs + 1):
    train_acc, train_loss = train(epoch)
    train_acc = test(train_loader, epoch=epoch)
    test_acc = test(test_loader, epoch, test=True)
    print('Epoch: {:03d}, Train Loss: {:.7f}, '
          'Train Acc: {:.7f}, Test Acc: {:.7f}'.format(epoch, train_loss,
                                                       train_acc, test_acc))
    
    y = np.zeros((len(test_dataset)))
    x = np.loadtxt(f"PROTEINS_epoch{epoch}.csv")
    for i in range(len(test_dataset)):
        y[i] = test_dataset[i].y
    y = torch.as_tensor(y).long()
    y_one_hot = F.one_hot(y, num_classes=dataset.num_classes).numpy()
    store_auc = 0
    for i in range(len(x[0,:])):
        auc = metrics.roc_auc_score(y_one_hot[:,i], x[:,i])
        print(f"AUC of {i} is: {auc}")
        store_auc += auc
    print(f"Average AUC: {store_auc/dataset.num_classes}")

print(f'Final Test Acc: {test_acc:.4f}')
