import os
import uproot
import numpy as np
import pandas as pd
import pickle as pkl

import torch
from torch.utils.data import TensorDataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from torchsummary import summary

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import RocCurveDisplay

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

seed = 7
np.random.seed(seed)

outputTag = '24-07-31'

outputNNDir =    'output/jetMatcherNNPyTorch/{}'.format(outputTag)
outputPlotsDir = 'plots/jetMatcherNNPyTorch/{}'.format(outputTag)
if not os.path.exists(outputNNDir): os.makedirs(outputNNDir)
if not os.path.exists(outputPlotsDir): os.makedirs(outputPlotsDir)

##############################
# LOAD & PREPARE DATA
##############################

batch_size = 50

print("reading datafile")
df = pd.read_csv('output/exportJetInfo/jets_uncompressed_le0p75mStop_24-07-30.csv')

masses = df.values[:,0:2]
inputs = df.values[:,2:-3]
nInputs = inputs.shape[1]
targets = df.values[:,-3:]

cutoff = masses[:, 1]/masses[:, 0] < 0.5
inputs = inputs[cutoff]
targets = targets[cutoff]

X_train, X_test, Y_train, Y_test = train_test_split(inputs,targets,test_size=0.1,random_state=seed, stratify=targets)
X_train, X_val, Y_train, Y_val   = train_test_split(X_train,Y_train,test_size=0.111,random_state=seed, stratify=Y_train)

scaler = StandardScaler().fit(X_train)
with open('{}/scaler.pkl'.format(outputNNDir),'wb') as f: pkl.dump(scaler,f)

X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
X_val = scaler.transform(X_val)

tXTrain = torch.Tensor(X_train)
tYTrain = torch.Tensor(Y_train)#.type(torch.LongTensor)
datasetTrain = TensorDataset(tXTrain,tYTrain)
loaderTrain = DataLoader(datasetTrain,batch_size=batch_size)

tXTest = torch.Tensor(X_test)
tYTest = torch.Tensor(Y_test)#.type(torch.LongTensor)
datasetTest = TensorDataset(tXTest,tYTest)
loaderTest = DataLoader(datasetTest,batch_size=100)

tXVal = torch.Tensor(X_val)
tYVal = torch.Tensor(Y_val)#.type(torch.LongTensor)
datasetVal = TensorDataset(tXVal,tYVal)
loaderVal = DataLoader(datasetVal,batch_size=100)

##############################
# DEFINE MODEL
##############################

learningRate = 0.001
print("defining net")
class Net(nn.Module):
  def __init__(self):
    super(Net,self).__init__()
    self.fc1 = nn.Linear(nInputs,nInputs)
    self.fc2 = nn.Linear(nInputs,3)
  def forward(self,x):
    x = F.relu(self.fc1(x))
    x = F.softmax(self.fc2(x),dim=1)
    return x

network = Net()
optimizer = optim.Adam(network.parameters(),lr=learningRate)

summary(network)

##############################
# DEFINE TRAIN & TEST
##############################

nEpochs = 100

#trainLosses = []
#trainCounter = []
#testLosses = []
#testCounter = [i*len(loaderTrain.dataset) for i in range(nEpochs + 1)]

lossFunction = nn.CrossEntropyLoss()

print("defining train, test, evaluate")
def train(epoch,loaderTrain):
  network.train()
  for iBatch,(data,target) in enumerate(loaderTrain):
    optimizer.zero_grad()
    output = network(data)
    loss = lossFunction(output,target)
    pred = output.data.max(1,keepdim=True)[1]
    target = target.data.max(1,keepdim=True)[1]
    nCorrect = pred.eq(target.data.view_as(pred)).sum()
    acc = nCorrect / len(data)
    loss.backward()
    optimizer.step()
    if iBatch % 20 == 0:
      print('Epoch {}: [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,iBatch * len(data),len(loaderTrain.dataset),100. * iBatch / len(loaderTrain),loss.item()))
      #trainLosses.append(loss.item())
      #trainCounter.append(iBatch * len(data) + ((epoch) * len(loaderTrain.dataset) / nEpochs))
  return loss.item(),acc.item()

def test(loaderTest):
  network.eval()
  lossTest = 0
  nCorrect = 0
  with torch.no_grad():
    for data,target in loaderTest:
      output = network(data)
      lossFunctionTest = nn.CrossEntropyLoss(reduction='sum')
      lossTest += lossFunctionTest(output,target).item()
      pred = output.data.max(1,keepdim=True)[1]
      target = target.data.max(1,keepdim=True)[1]
      nCorrect += pred.eq(target.data.view_as(pred)).sum()
  lossTest /= len(loaderTest.dataset)
  acc = nCorrect / len(loaderTest.dataset)
  #testLosses.append(lossTest)
  print('Test: Average loss = {:.4f}, Accuracy = {}/{} ({:.0f}%)'.format(lossTest,nCorrect,len(loaderTest.dataset),100. * acc))
  return lossTest,acc

def evaluate(loaderVal):
  network.eval()
  lossTest = 0
  nCorrect = 0
  outputs = []
  targets = []
  with torch.no_grad():
    for data,target in loaderTest:
      targets.extend(target.data.tolist())
      output = network(data)
      lossFunctionTest = nn.CrossEntropyLoss(reduction='sum')
      lossTest += lossFunctionTest(output,target).item()
      pred = output.data.max(1,keepdim=True)[1]
      target = target.data.max(1,keepdim=True)[1]
      nCorrect += pred.eq(target.data.view_as(pred)).sum()
      outputs.extend(output.data.tolist())
  lossTest /= len(loaderTest.dataset)
  acc = nCorrect / len(loaderTest.dataset)
  #testLosses.append(lossTest)
  print('Test: Average loss = {:.4f}, Accuracy = {}/{} ({:.0f}%)'.format(lossTest,nCorrect,len(loaderTest.dataset),100. * acc))
  return lossTest,acc,outputs,targets

##############################
# TRAIN
##############################
print("training")
lossesTrain = []
accTrain = []
lossesTest = []
accTest = []

for epoch in range(nEpochs):
  loss,acc = train(epoch,loaderTrain)
  lossesTrain.append(loss)
  accTrain.append(acc)
  loss,acc = test(loaderTest)
  lossesTest.append(loss)
  accTest.append(acc)
print('Training complete. Running validation...')
loss,acc,outputs,targets = evaluate(loaderVal)
#print(trainCounter)
#print(lossesTrain)
#print(testCounter)
#print(lossesTest)

torch.save(network,'{}/jetMatcherNN.pt'.format(outputNNDir))

plt.figure(figsize=(15,10))

rocLabels = ['Stop vs. Rest','Chargino vs. Rest','Other vs. Rest']
for i in [0,1,2]:
  fpr, tpr, thresholds = roc_curve(np.array(targets)[:,i],np.array(outputs)[:,i])
  roc_auc = auc(fpr, tpr)
  rocDisplay = RocCurveDisplay(fpr=fpr,tpr=tpr,roc_auc=roc_auc,estimator_name=rocLabels[i])
  #RocCurveDisplay.from_predictions(Y_test_onehot[:,0],Y_pred,name='stop vs. rest',color='red')
  ax = plt.subplot(2,2,3)
  rocDisplay.plot(ax=ax)
ax.plot([0,1],[0,1],'k--',label='chance (AUC = 0.5)')
ax.set(xlabel='False Positive Rate',ylabel='True Positive Rate')
ax.legend()

confusionMatrix = confusion_matrix(np.argmax(targets,1),np.argmax(outputs,1))
cmPlot = ConfusionMatrixDisplay(confusion_matrix=confusionMatrix,display_labels=['stop','chargino','other'])
cmPlot.plot(ax=plt.subplot(2,2,4))

ax = plt.subplot(2,2,1)
ax.plot(np.linspace(1,nEpochs,nEpochs),lossesTrain,label='Training Loss')
ax.plot(np.linspace(1,nEpochs,nEpochs),lossesTest,label='Test Loss')
ax.legend(loc='upper right')
ax.set_xlabel('Epoch')
ax.set_ylabel('Loss')

ax = plt.subplot(2,2,2)
ax.plot(np.linspace(1,nEpochs,nEpochs),accTrain,label='Training Accuracy')
ax.plot(np.linspace(1,nEpochs,nEpochs),accTest,label='Test Accuracy')
ax.legend(loc='lower right')
ax.set_xlabel('Epoch')
ax.set_ylabel('Accuracy')

plt.savefig('{}/performance.pdf'.format(outputPlotsDir))
