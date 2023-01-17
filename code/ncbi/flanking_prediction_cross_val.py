#import dependencies
import numpy as np
import pandas as pd
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report,roc_curve,auc,RocCurveDisplay

#import features
features=pd.read_csv('../data/flanking_features.csv',index_col=0)

#create labels and add them to the dataframe
labels=pd.Series(features.index).apply(lambda name: 0 if name[0]=='r' else 1)
features['label']=labels.values

#separate labels
X,y=features.drop('label',axis=1),features['label']

cv = StratifiedKFold(n_splits=10)
classifier=RandomForestClassifier(random_state=1,n_estimators=20,max_depth=5)

#true positives and aucss
tprs,aucs = [],[]
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots()
for i, (train, test) in enumerate(cv.split(X, y)):
    classifier.fit(X.iloc[train], y.iloc[train])
    viz = RocCurveDisplay.from_estimator(
        classifier,
        X.iloc[test],
        y.iloc[test],
        name=f"ROC fold {i}",
        alpha=0.2,
        lw=1,
        ax=ax,
    )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    '-',
    color="b",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=1,
)

#visualize std
std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(
    mean_fpr,
    tprs_lower,
    tprs_upper,
    color="grey",
    alpha=0.2,
    label=r"$\pm$ 1 std. dev.",
)

ax.set(
    xlim=[-0.05, 1.05],
    ylim=[-0.05, 1.05],
    xlabel='TPR',
    ylabel='FPR'
)
ax.plot([0,1],[0,1],'--',lw=2,label='random')

ax.set_xlabel('FPR',fontsize=20),
ax.set_ylabel('TPR',fontsize=20)

ax.legend(loc="lower right")

plt.tight_layout()

plt.savefig('../results/flanking_cvrocs.png',dpi=400)

#export feature importances
feature_importances=pd.Series(classifier.feature_importances_)
feature_importances.index=X.columns
feature_importances.to_csv('../results/flanking_feature_importances.csv')