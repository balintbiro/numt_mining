#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from dask import dataframe as dd
from dask_ml.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from dask_ml.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

#read features
features=dd.read_csv('../data/flanking_features.csv')

#create labels and add them to the dataframe
labels=features['Unnamed: 0'].apply(lambda name: 0 if name[0]=='r' else 1)
features['label']=labels

#update features; delete fasta headers (genomic ids)
features=features.drop('Unnamed: 0',axis=1)

#separate labels and data
X,y=features.drop('label',axis=1),features['label']

#scale X
scaler=StandardScaler()
X_scaled=scaler.fit_transform(X)

#create training and test sets
X_train,X_test,y_train,y_test=train_test_split(X_scaled,y,random_state=0,shuffle=False)

#initialize classifier
clf=RandomForestClassifier(random_state=0)
clf.fit(X_train,y_train)
print(clf.score(X_test,y_test))