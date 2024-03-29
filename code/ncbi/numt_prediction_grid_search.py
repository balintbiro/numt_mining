#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split,RandomizedSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

#import features
features=pd.read_csv('../data/iFeatureOmegaCLI_features.csv',index_col=0)
features=features.sample(n=10000,replace=False)

#separate labels
X,y=features.drop(['label','order','order_label'],axis=1),features['label']

#scale X
X_scaled=StandardScaler().fit_transform(X)

#split training and test set
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, random_state=0)

#initialize model
rfc = RandomForestClassifier(random_state=0)

#setting parameters
param_dist = {
    'max_depth': [1,2,3],#max depth of a tree
    'max_features': [1,2,3,4,10,25,50],#The number of features to consider when looking for the best split:
    'min_samples_leaf': np.linspace(1,10,5,dtype=int),#The minimum number of samples required to be at a leaf node
    'min_samples_split': np.linspace(1,30,5,dtype=int),#The minimum number of samples required to split an internal node
    'n_estimators': np.linspace(1,10,5,dtype=int)#number of trees
}

#setting grid search for hyperparameter optimisation
grid_search = RandomizedSearchCV(estimator = rfc, param_distributions = param_dist,n_jobs=-1,
                          verbose = 0,cv=10,n_iter=100,scoring='roc_auc')

#grid search for hyperparameter optimisation
grid_search.fit(X_train, y_train)

#transform the reults and save them
gsCV_results=pd.DataFrame.from_dict(grid_search.cv_results_)
gsCV_results.to_csv('../results/gsCV_results.csv')

####################################################################################################
#                                       visualisation                                              #
#           for just the visualization, copy the code from here plus the dependencies              #
####################################################################################################

#read in results
gscv_res=pd.read_csv('../results/gsCV_results.csv',index_col=0)

#filter columns
params=list(gscv_res.filter(like='param_'))
params.append('mean_test_score')
gscv_res=rgscv_res[params]
param_dict=pd.Series({
  'param_n_estimators':'n estimators',
  'param_min_samples_split':'mn split',
  'param_min_samples_leaf':'mn leaf',
  'param_max_features':'mx features',
  'param_max_depth':'mx depth',
  'mean_test_score':'avg auc'
})
gscv_res.columns=param_dict[params].values

#create plot
fig = px.parallel_coordinates(gscv_res, color="avg auc",
                              dimensions=gscv_res.columns,
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
    '../results/grid_search.png',
    validate=True,
    width=650,
    height=350,
    scale=3
  )