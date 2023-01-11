#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

#import features
features=pd.read_csv('../data/iFeatureOmegaCLI_features.csv',index_col=0)

#separate labels
X,y=features.drop(['label','order','order_label'],axis=1),features['label']

#scale X
X_scaled=StandardScaler().fit_transform(X)

#split training and test set
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, random_state=0)

#initialize model
rfc = RandomForestClassifier(random_state=0)

#setting parameters
param_grid = {
    'max_depth': [2,3,4,5],#max depth of a tree
    'max_features': np.linspace(2,len(X.columns),5,dtype=int),#The number of features to consider when looking for the best split:
    'min_samples_leaf': np.linspace(2,10,5,dtype=int),#The minimum number of samples required to be at a leaf node
    'min_samples_split': np.linspace(2,100,5,dtype=int),#The minimum number of samples required to split an internal node
    'n_estimators': np.linspace(3,20,5,dtype=int)#number of trees
}

#setting grid search for hyperparameter optimisation
grid_search = GridSearchCV(estimator = rfc, param_grid = param_grid,n_jobs=-1,
                          verbose = 0,scoring='roc_auc')

#grid search for hyperparameter optimisation
grid_search.fit(X_train, y_train)

#transform the reults and save them
gsCV_results=pd.DataFrame.from_dict(grid_search.cv_results_)
gsCV_results.to_csv('../results/gsCV_results.csv')