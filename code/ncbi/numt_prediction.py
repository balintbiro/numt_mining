#import dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve

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
       'downstream_entropy', 'label']]

#separate labels
X=features.loc[:,features.columns!='label']
y=features['label']

#split training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

#initialize model
model = RandomForestClassifier()

#train model
model.fit(X_train,y_train)

#prediction on test set
y_pred_test=model.predict(X_test)

print(f'Accuracy score is {accuracy_score(y_test, y_pred_test)}!')

#dataprep for ROC curve
fpr, tpr, _ = roc_curve(y_test,  y_pred_test)
fig,axs=plt.subplots(1,1)
axs.plot(fpr,tpr)
plt.savefig('../results/numt_rf_roc.png',dpi=200)
