#commonly cnn-s are used for image related tasks (classification, object detection, segmentation etc.)
import os
import cv2
import numpy as np
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

REBUILD_DATA=False#set ti True to once, then leave it False

class DogsVSCats():
	IMG_SIZE=50#thats gonna resize our images
	CATS="../data/kagglecatsanddogs_5340/PetImages/Cat/"
	DOGS="../data/kagglecatsanddogs_5340/PetImages/Dog/"
	LABELS={CATS:0,DOGS:0}
	training_data=[]
	catcount,dogcount=0,0
	
	def make_training_data(self):
		for label in self.LABELS:
			print(label)
			for f in tqdm(os.listdir(label)):
				try:
					path=os.path.join(label,f)
					img=cv2.imread(path,cv2.IMREAD_GRAYSCALE)#grayscale conversion
					img=cv2.resize(img,(self.IMG_SIZE,self.IMG_SIZE))
					self.training_data.append([np.array(img),np.eye(2)[self.LABELS[label]]])#one hot encoding

					if label==self.CATS:
						self.catcount+=1
					elif label==self.DOGS:
						self.dogcount+=1
				except Exception as e:
					print(str(e))
					pass
		np.random.shuffle(self.training_data)
		np.save("../data/training_data.npy",self.training_data)
		print("Cats: ",self.catcount)
		print("Dogs: ",self.dogcount)


if REBUILD_DATA:
	dogsvcats=DogsVSCats()
	dogsvcats.make_training_data()

training_data=np.load("../data/training_data.npy",allow_pickle=True)

class Net(nn.Module):
	def __init__(self):
		super().__init__()
		self.conv1=nn.Conv2d(1,32,5)#last number is the kernel size
		self.conv2=nn.Conv2d(32,64,5)
		self.conv3=nn.Conv2d(64,128,5)
		

		x=torch.randn(50,50).view(-1,1,50,50)
		self._to_linear=None
		self.convs(x)
		self.fc1=nn.Linear(self._to_linear,512)
		self.fc2=nn.Linear(512,2)
		
	def convs(self,x):
		x=F.max_pool2d(F.relu(self.conv1(x)),(2,2))
		x=F.max_pool2d(F.relu(self.conv2(x)),(2,2))
		x=F.max_pool2d(F.relu(self.conv3(x)),(2,2))
		print(x[0].shape)

		if self._to_linear is None:
			self._to_linear=x[0].shape[0]*x[0].shape[1]*x[0].shape[2]
		return x

	def forward(self,x):
		x=self.convs(x)
		x=x.view(-1,self._to_linear)
		x=F.relu(self.fc1(x))
		x=self.fc2(x)
		return F.softmax(x,dim=1)

net=Net()
optimizer=optim.Adam(net.parameters(),lr=.001)
loss_function=nn.MSELoss()

X=torch.Tensor([i[0] for i in training_data]).view(-1,50,50)
X=X/255.0
y=torch.Tensor([i[1] for i in training_data])

VAL_PCT=.1
val_size=int(len(X)*VAL_PCT)
X_train,X_test=X[:-val_size],X[-val_size:]
y_train,y_test=y[:-val_size],y[-val_size:]

BATCH_SIZE=100
EPOCHS=1

for epoch in range(EPOCHS):
	for i in tqdm(range(0,len(X_train),BATCH_SIZE)):
		batch_X=X_train[i:i+BATCH_SIZE].view(-1,1,50,50)
		batch_y=y_train[i:i+BATCH_SIZE]

		net.zero_grad()
		outputs=net(batch_X)
		
		loss=loss_function(outputs,batch_y)
		loss.backward()
		optimizer.step()

correct,total=0,0

with torch.no_grad():
	for i in tqdm(range(len(X_test))):
		real_class=torch.argmax(y_test[i])
		net_out=net(X_test[i].view(-1,1,50,50))[0]
		predicted_class=torch.argmax(net_out)
		if predicted_class==real_class:
			correct+=1
		total+=1

print('Accuracy: ',round(correct/total,3))
