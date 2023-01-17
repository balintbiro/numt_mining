#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split, RandomizedSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

#read features
features=pd.read_csv('../data/flanking_features.csv',index_col=0)
features=features.sample(n=10000,replace=False)

#create labels and add them to the dataframe
labels=pd.Series(features.index).apply(lambda name: 0 if name[0]=='r' else 1)
features['label']=labels.values

#separate labels and data
X,y=features.drop('label',axis=1),features['label']

#scale X
scaler=StandardScaler()
X_scaled=scaler.fit_transform(X)

#create training and test sets
X_train,X_test,y_train,y_test=train_test_split(X_scaled,y,random_state=0)

#initialize classifier
rfc=RandomForestClassifier(random_state=0)

#setting parameters
param_dist = {
    'max_depth': [2,3,4,5],#max depth of a tree
    'max_features': np.linspace(2,100,5,dtype=int),#The number of features to consider when looking for the best split:
    'min_samples_leaf': np.linspace(2,10,5,dtype=int),#The minimum number of samples required to be at a leaf node
    'min_samples_split': np.linspace(2,30,5,dtype=int),#The minimum number of samples required to split an internal node
    'n_estimators': np.linspace(3,10,5,dtype=int)#number of trees
}

#setting grid search for hyperparameter optimisation
randomCV = RandomizedSearchCV(estimator = rfc, param_distributions = param_dist,n_jobs=-1,
                          verbose = 0,n_iter=100,scoring='roc_auc')

#grid search for hyperparameter optimisation
randomCV.fit(X_train, y_train)

#transform the reults and save them
randomCV_results=pd.DataFrame.from_dict(randomCV.cv_results_)
randomCV_results.to_csv('../results/flanking_gsCV_results.csv')

####################################################################################################
#                                       visualisation                                              #
#           for just the visualization, copy the code from here plus the dependencies              #
####################################################################################################

#read in results
randomcv_res=pd.read_csv('../results/flanking_gsCV_results.csv',index_col=0)

#filter columns
params=list(randomcv_res.filter(like='param_'))
params.append('mean_test_score')
randomcv_res=randomcv_res[params]
randomcv_res.columns=['mx depth','mx features','mn leaf','mn split','n estimators','avg auc']

#create plot
fig = px.parallel_coordinates(gscv_res, color="avg auc",
                              dimensions=randomcv_res.columns,
                              color_continuous_scale=px.colors.diverging.Tealrose,
                              color_continuous_midpoint=.5)

#update layout
fig.update_layout(
    height=350,
        width=650,
    font=dict(
        family='Arial',
        size=20,
        color='#000000'
      )
  )

#create figure
fig.write_image(
    '../results/flanking_grid_search.png',
    validate=True,
    width=650,
    height=350,
    scale=3
  )