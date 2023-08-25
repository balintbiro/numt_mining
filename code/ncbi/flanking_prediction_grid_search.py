#import dependencies
import os
import xgboost
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
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
classifier=xgboost.XGBClassifier()

#setting parameters
param_dist = {
    'learning_rate':[.05,.1,.15,.2,.25,.3],
    'max_depth':[3,4,5,6,8,10,12,15],
    'min_child_weight':[1,3,5,7],
    'gamma':[0,.1,.2,.3,.4],
    'colsample_bytree':[.3,.4,.5,.7]
}

#setting grid search for hyperparameter optimisation
randomCV = RandomizedSearchCV(estimator = classifier, param_distributions = param_dist,n_jobs=-1,
                          verbose = 0,cv=10,n_iter=10,scoring='roc_auc')

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
param_dict=pd.Series({
  'param_n_estimators':'n estimators',
  'param_min_samples_split':'mn split',
  'param_min_samples_leaf':'mn leaf',
  'param_max_features':'mx features',
  'param_max_depth':'mx depth',
  'mean_test_score':'avg auc'
})
randomcv_res.columns=param_dict[params].values

#create plot
fig = px.parallel_coordinates(randomcv_res, color="avg auc",
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