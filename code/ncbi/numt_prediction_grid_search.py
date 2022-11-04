#import dependencies
import numpy as np
import pandas as pd
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

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

#avoid imbalance-sample randomly
features=pd.concat([
    features[features['label']==1].sample(n=len(features[features['label']==0]),replace=False),
    features[features['label']==0]
    ])

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
param_grid = {
    'max_depth': np.linspace(1,100,10,dtype=int),#max depth of a tree
    'max_features': np.linspace(1,len(X.columns),10,dtype=int),#The number of features to consider when looking for the best split:
    'min_samples_leaf': np.linspace(1,10,3,dtype=int),#The minimum number of samples required to be at a leaf node
    'min_samples_split': np.linspace(1,100,10,dtype=int),#The minimum number of samples required to split an internal node
    'n_estimators': np.linspace(10,1000,10,dtype=int)#number of trees
}

#setting grid search for hyperparameter optimisation
grid_search = GridSearchCV(estimator = rfc, param_grid = param_grid, 
                          cv = 5, n_jobs = -1, verbose = 2,scoring='roc_auc')

#grid search for hyperparameter optimisation
grid_search.fit(X_train, y_train)

#transform the reults and save them
gsCV_results=pd.DataFrame.from_dict(grid_search.cv_results_)

#select parameters and test score
fil=pd.Series(gsCV_results).str.contains('param_m|param_n|mean_test_score').tolist()

#visualize result
fig,axs=plt.subplots(1,1)
plt.style.use('ggplot')
sns.scatterplot(x='param_max_depth',
               y='mean_test_score',
               hue='param_max_features',
                size='param_min_samples_leaf',
                style='param_n_estimators',
                palette='tab10',
               data=cv_res.loc[:,fil].dropna(),
               ax=axs)
plt.xlabel('Max depth (log)',fontsize=20)
plt.ylabel('Test AUC',fontsize=20)
plt.xscale('log')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('../results/gs_results.png',transparent=False,dpi=200)