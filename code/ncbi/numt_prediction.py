#import dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import StratifiedKFold
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc

#import features
features=pd.read_csv('../data/ml_features.csv')

#clear features df
features=features.dropna()

#filter df
fil=features.apply(lambda row:(row['upstream_size']>4900) and (row['downstream_size']>4900),axis=1)
features=features[fil]
features=features[['GC', 'upstream_GC', 'downstream_GC', 'uSW_mean', 'uSW_median',
       'uRMs_count', 'uRMs_lengths', 'dSW_mean', 'dSW_median', 'dRMs_count',
       'dRMs_lengths', 'rel_start', 'entropy', 'upstream_entropy',
       'downstream_entropy', 'label','u_TmGC', 'u_TmNN', 'u_TmW', 'TmGC', 'TmNN', 'TmW', 'd_TmGC', 'd_TmNN',
       'd_TmW']]

#just the sequence itself
#features=features[['GC', 'rel_start', 'entropy', 'label','TmGC', 'TmNN', 'TmW',]]

#separate labels
X=features.loc[:,features.columns!='label']
y=features['label']

#split training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

#initialize model
rfc = RandomForestClassifier()

#train model
rfc.fit(X_train,y_train)

#setting parameters
#param_grid = {
#    'bootstrap': [True,False],
#    'max_depth': np.linspace(1,10,5,dtype=int),
#    'max_features': [1,3,5],
#    'min_samples_leaf': np.linspace(1,10,3,dtype=int),
#    'min_samples_split': np.linspace(1,100,3,dtype=int),
#    'n_estimators': np.linspace(10,1000,3,dtype=int)
#}

#setting grid search for hyperparameter optimisation
#grid_search = GridSearchCV(estimator = rfc, param_grid = param_grid, 
#                          cv = 5, n_jobs = -1, verbose = 2,scoring='roc_auc')

#grid search for hyperparameter optimisation
#grid_search.fit(X_train, y_train)

#transform the reults and save them
#gsCV_results=pd.DataFrame.from_dict(grid_search.cv_results_)
#gsCV_results.to_csv('../results/gsCV_results.csv')