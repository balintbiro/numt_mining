import numpy as np
import pandas as pd
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay,roc_auc_score

import xgboost
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_validate

df=pd.read_csv(
    '../data/iFeatureOmegaCLI_features.csv',
    index_col=0
)

X,y=df.drop(columns=['label','order_label','order']),df.label

def cross_validation(classifier_tuple):
	global cv_results
	results=cross_validate(classifier_tuple[0],X,y,cv=10,scoring=('accuracy','roc_auc'))
	results=pd.DataFrame(results)
	results['model']=classifier_tuple[1]
	cv_results.append(results)

classifiers=pd.Series(
    data=[
    	#(SVC(),'SVC'),
        (xgboost.XGBClassifier(),'XGB'),
        #(RandomForestClassifier(),'RFC'),
        #(LogisticRegression(),'LR'),
        #(MLPClassifier(),'MLP'),
        #(GaussianNB(),'GNB')
    ]
)

cv_results=[]
classifiers.apply(cross_validation)

cv_results=pd.concat(cv_results)

cv_results.to_csv('../data/ml_model_testing_results.csv')


####################################################################################################
#                                       visualisation                                              #
#           for just the visualization, copy the code from here plus the dependencies              #
####################################################################################################

df=pd.read_csv('../../results/ml_model_testing_results.csv',index_col=0)

grouped_df=df.groupby(by='model').agg(['mean','std'])
grouped_df.columns=pd.Series(grouped_df.columns.values).str.join(sep='_').values

fig,axs=plt.subplots()
x,y=grouped_df['test_roc_auc_mean'].values,grouped_df['test_accuracy_mean'].values
model_scatter=axs.scatter(x,y)
axs.errorbar(
    x=grouped_df['test_roc_auc_mean'].values,y=grouped_df['test_accuracy_mean'].values,
    xerr=grouped_df['test_roc_auc_std'].values,yerr=grouped_df['test_accuracy_std'].values,
    fmt='o',c='g'
)
axs.annotate(grouped_df.index[0],xy=(x[0]+.02,y[0]-0.01),textcoords='offset points',fontsize=15)
axs.annotate(grouped_df.index[1],xy=(x[1]-.03,y[1]-0.01),textcoords='offset points',fontsize=15)
axs.annotate(grouped_df.index[2],xy=(x[2]-.03,y[2]-0.01),textcoords='offset points',fontsize=15)
axs.annotate(grouped_df.index[3],xy=(x[3]-.03,y[3]-0.01),textcoords='offset points',fontsize=15)
axs.annotate(grouped_df.index[4],xy=(x[4]-.03,y[4]-0.01),textcoords='offset points',fontsize=15)
axs.annotate(grouped_df.index[5],xy=(x[5]-.03,y[5]-0.01),textcoords='offset points',fontsize=15)

axs.set_ylabel('Accuracy',fontsize=20)
axs.set_xlabel('AUC',fontsize=20)
plt.tight_layout()
plt.savefig('../../results/ml_testing.png',dpi=400)