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
features=pd.read_csv('../data/iFeatureOmegaCLI_features.csv',index_col=0)

#separate labels
X,y=features.drop(['label','order','order_label'],axis=1),features['label']

cv = StratifiedKFold(n_splits=10)
classifier=RandomForestClassifier(random_state=1,n_estimators=7,min_samples_leaf=1,max_features=50,max_depth=3)

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

plt.savefig('../results/cvrocs.png',dpi=400)

#export feature importances
feature_importances=pd.Series(classifier.feature_importances_)
feature_importances.index=X.columns
feature_importances.to_csv('../results/feature_importances.csv')

#visualize feature importances
####################################################################################################
#                                       visualisation                                              #
#           for just the visualization, copy the code from here plus the dependencies              #
####################################################################################################

#import data
feat_imp=pd.read_csv('../results/feature_importances.csv',index_col=0)
feat_imp.columns=['importance']
feat_imp['desc_group']=(64*['Kmer Type1'])+(32*['RCKmer Type1'])+(32*['RCKmer Type2'])+(64*['Mismatch'])+(9*['Zcurve'])+(18*['NMBroto'])+(64*['Kmer Type2'])+(4*['NAC'])

#sort every descriptor group
feat_groups=feat_imp['desc_group'].unique()
importances=pd.Series(feat_groups).apply(lambda feat_group: np.array(feat_imp.loc[feat_imp['desc_group']==feat_group]['importance']))
importances.index=feat_groups
importances.apply(lambda imp_lst: imp_lst.sort())

#define colors
colors=['blue','grey','red','pink','green','black','orange','brown',]

#plotter function
def plotter(imp_array):
    global x,indexer
    axs.barh(np.arange(x,(x+len(imp_array))),imp_array,height=1,color=colors[indexer])
    axs.fill_between((0,0.1),x,(x+len(imp_array)),color=colors[indexer],alpha=.05)
    axs.text(0.065,(x+len(imp_array)/2),importances.index[indexer],fontsize=15)
    x+=len(imp_array)
    indexer+=1

fig,axs=plt.subplots(1,1,figsize=(6,10))
axs.set_xlim(0,.1)
x,indexer=0,0
importances.apply(plotter)
axs.set(yticklabels=[],yticks=[])
axs.set_ylabel('Features',fontsize=20)
axs.set_xlabel('Feature importance',fontsize=20)
plt.tight_layout()
plt.savefig('../results/feature_importances.png',dpi=400)