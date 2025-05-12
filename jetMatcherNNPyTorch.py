import os
import uproot
import numpy as np
import pandas as pd
import pickle as pkl
from datetime import datetime
import argparse

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

timeTag = datetime.now().strftime("%y-%m-%d-%H-%M")

seed = 7
np.random.seed(seed)

parser = argparse.ArgumentParser(description='Jet Matcher NN PyTorch')
parser.add_argument('--input',type=str,required=True,help='Path to input file')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('--test',action='store_true',default=False,help='Do not save outputs')
args = parser.parse_args()

# Input
if not os.path.exists(args.input): 
  print('ERROR: Provided input file does not exist. Exiting...')
  quit()

# Output
if not args.test:
  outputDir = 'output/jetMatcherNNPyTorch/{}_{}'.format(args.tag,timeTag)
  if not os.path.exists(outputDir): os.makedirs(outputDir)

##############################
# LOAD & PREPARE DATA
##############################

batch_size = 50

df = pd.read_csv(args.input)

# Select training points
print('{} total jets for filtering'.format(len(df.index)))
df = df[(df['mChargino'] / df['mStop']) > (2/3)]
print('{} total jets after filtering'.format(len(df.index)))

# Inputs
#inputs = df.drop(['mStop','mChargino','isStopMatched','isChiMatched','isOther'],axis=1).values
inputs = df[['jetOrdinality','jetPT','jetEta','jetPhi','jetBScore','m3M','m3PT','m3Eta','m3Phi','m4M','m4PT','m4Eta','m4Phi','nJets']].values
#inputs = df[['jetOrdinality','jetPT','jetEta','jetPhi','jetBScore','nJets']].values
nInputs = inputs.shape[1]

# Outputs
#targets = df[['isStopMatched','isChiMatched','isOther']].values
targets = df['isChiMatched']
nOutputs = 1

# Number of nodes in hidden layer = average of input and output nodes
nNodesHidden = np.floor((nInputs + nOutputs)/2).astype(int)

# Set weights to balance different signal points
df['mStop'] = df['mStop'].astype(str)
df['mChargino'] = df['mChargino'].astype(str)
df['mStop_mChargino'] = df[['mStop','mChargino']].agg('_'.join,axis=1)
weightVals = 1 / df['mStop_mChargino'].value_counts(normalize=True)
weights = df['mStop_mChargino'].map(weightVals)

X_train, X_test, Y_train, Y_test, W_train, W_test   = train_test_split(inputs,targets,weights,test_size=0.1,random_state=seed,stratify=targets)
X_train, X_val, Y_train, Y_val, W_train, W_val      = train_test_split(X_train,Y_train,W_train,test_size=0.111,random_state=seed,stratify=Y_train)

scaler = StandardScaler().fit(X_train)
if not args.test: 
  with open('{}/scaler.pkl'.format(outputDir),'wb') as f: pkl.dump(scaler,f)

X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
X_val = scaler.transform(X_val)

tXTrain = torch.Tensor(X_train)
tYTrain = torch.Tensor(Y_train.to_numpy())
tWTrain = torch.Tensor(W_train.to_numpy())
datasetTrain = TensorDataset(tXTrain,tYTrain,tWTrain)
loaderTrain = DataLoader(datasetTrain,batch_size=batch_size)

tXTest = torch.Tensor(X_test)
tYTest = torch.Tensor(Y_test.to_numpy())
tWTest = torch.Tensor(W_test.to_numpy())
datasetTest = TensorDataset(tXTest,tYTest,tWTest)
loaderTest = DataLoader(datasetTest,batch_size=100)

tXVal = torch.Tensor(X_val)
tYVal = torch.Tensor(Y_val.to_numpy())
tWVal = torch.Tensor(W_val.to_numpy())
datasetVal = TensorDataset(tXVal,tYVal,tWVal)
loaderVal = DataLoader(datasetVal,batch_size=100)

##############################
# DEFINE MODEL
##############################

learningRate = 0.001 # 0.00005 Lower rate

class Net(nn.Module):
  def __init__(self):
    super(Net,self).__init__()
    self.fc1 = nn.Linear(nInputs,nNodesHidden)
    self.fc2 = nn.Linear(nNodesHidden,1)
  def forward(self,x):
    x = F.relu(self.fc1(x))
    x = F.sigmoid(self.fc2(x)) #F.softmax(self.fc2(x),dim=1)
    return x

torch.manual_seed(seed)
network = Net()
optimizer = optim.Adam(network.parameters(),lr=learningRate)

summary(network)

##############################
# DEFINE TRAIN & TEST
##############################

nEpochs = 50 # 150

#lossFunction = nn.CrossEntropyLoss(reduction='none')
lossFunction = nn.BCELoss(reduction='none')

def weightedLoss(output,target,weight):
  weight = 1
  return (lossFunction(output,target) * weight).mean()

def train(epoch,loaderTrain):
  network.train()
  lossSum = 0
  nCorrect = 0
  for iBatch,(data,target,weights) in enumerate(loaderTrain):
    optimizer.zero_grad()
    output = torch.squeeze(network(data))
    #loss = lossFunction(output,target)
    loss = weightedLoss(output,target,weights)
    lossSum += loss
    pred = torch.round(output) #output.data.max(1,keepdim=True)[1]
    target = target #target.data.max(1,keepdim=True)[1]
    nCorrect += pred.eq(target.data.view_as(pred)).sum()
    loss.backward()
    optimizer.step()
    if iBatch % 200 == 0:
      print('Epoch {}: [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch+1,iBatch * len(data),len(loaderTrain.dataset),100. * iBatch / len(loaderTrain),loss.item()))
  lossAvg = lossSum / len(loaderTrain)
  acc = nCorrect / len(loaderTrain.dataset)
  return lossAvg.item(),acc.item()

def test(loaderTest):
  network.eval()
  loss = 0
  nCorrect = 0
  with torch.no_grad():
    for (data,target,weights) in loaderTest:
      output = torch.squeeze(network(data))
      #lossFunctionTest = nn.CrossEntropyLoss(reduction='sum')
      #loss += lossFunctionTest(output,target).item()
      loss += weightedLoss(output,target,weights)
      pred = torch.round(output) #output.data.max(1,keepdim=True)[1]
      target = target #target.data.max(1,keepdim=True)[1]
      nCorrect += pred.eq(target.data.view_as(pred)).sum()
  loss /= len(loaderTest)
  acc = nCorrect / len(loaderTest.dataset)
  return loss,acc

def evaluate(loaderVal):
  network.eval()
  lossTest = 0
  nCorrect = 0
  outputs = []
  targets = []
  with torch.no_grad():
    for (data,target,weights) in loaderTest:
      targets.extend(target.data.tolist())
      output = torch.squeeze(network(data))
      #lossFunctionTest = nn.CrossEntropyLoss(reduction='sum')
      #lossTest += lossFunctionTest(output,target).item()
      lossTest += weightedLoss(output,target,weights)
      pred = torch.round(output) #output.data.max(1,keepdim=True)[1]
      target = target #target.data.max(1,keepdim=True)[1]
      nCorrect += pred.eq(target.data.view_as(pred)).sum()
      outputs.extend(output.data.tolist())
  lossTest /= len(loaderTest)
  acc = nCorrect / len(loaderTest.dataset)
  return lossTest,acc,outputs,targets

##############################
# TRAIN
##############################

lossesTrain = []
accTrain = []
lossesTest = []
accTest = []

for epoch in range(nEpochs):
  loss,acc = train(epoch,loaderTrain)
  lossesTrain.append(loss)
  accTrain.append(acc)
  print('----- Epoch {} Train/Test Stats -----'.format(epoch+1))
  print('Train:\tAverage loss = {:.5f},\tAccuracy = {:.1f}%'.format(loss,acc * 100.))
  loss,acc = test(loaderTest)
  lossesTest.append(loss)
  accTest.append(acc)
  print('Test:\tAverage loss = {:.5f},\tAccuracy = {:.1f}%'.format(loss,acc * 100.))
  print('-------------------------------------{}'.format('-' if epoch >= 10 else ''))
print('Training complete. Running validation...')
loss,acc,outputs,targets = evaluate(loaderVal)

if not args.test: 
  print('Saving model to {}'.format(outputDir))
  torch.save(network,'{}/jetMatcherNN.pt'.format(outputDir))

plt.figure(figsize=(15,10))

rocLabels = ['Stop vs. Rest','Chargino vs. Rest','Other vs. Rest']
isBinary = isinstance(targets[0],float)
for i in [1]: #[0,1,2]:
  fpr, tpr, thresholds = roc_curve( np.array(targets)[:,i] if not isBinary else np.array(targets),
                                    np.array(outputs)[:,i] if not isBinary else np.array(outputs))
  roc_auc = auc(fpr, tpr)
  print('{} AUC = {}'.format(rocLabels[i],roc_auc))
  rocDisplay = RocCurveDisplay(fpr=fpr,tpr=tpr,roc_auc=roc_auc,estimator_name=rocLabels[i])
  #RocCurveDisplay.from_predictions(Y_test_onehot[:,0],Y_pred,name='stop vs. rest',color='red')
  ax = plt.subplot(2,2,3)
  rocDisplay.plot(ax=ax)
ax.plot([0,1],[0,1],'k--',label='chance (AUC = 0.5)')
ax.set(xlabel='False Positive Rate',ylabel='True Positive Rate')
ax.legend()

if not isBinary:
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

if not args.test: plt.savefig('{}/performance.pdf'.format(outputDir))
