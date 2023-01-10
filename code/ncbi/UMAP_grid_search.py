#import dependencies
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

#read in dataframe
features=pd.read_csv('../data/iFeatureOmegaCLI_features.csv',index_col=0)

#add label names
label_dict=pd.Series(['NUMT','random'])
label_dict.index=[1,0]
features['label_name']=label_dict[features['label']].values

#for order UMAP
features=features.loc[features['order_label']!=-1].dropna()

#sample df
#features=features.sample(n=int(len(features)/40),replace=False)

X=features.drop(['label','order','order_label','label_name'],axis=1)

#define function for plotting the result
def plotter(curr_ax,title=None):
    kwargs={'edgecolor':'face'}
    scplot=sns.scatterplot(x='x',y='y',hue='order',data=features,s=10,palette='tab20',alpha=0.25,**kwargs)
    curr_ax.set(xticklabels=[],yticklabels=[])
    #scplot.legend(title='Sequence',fontsize=15,title_fontsize=15)#this is for the numt vs random umap
    curr_ax.set_xlabel('UMAP1',fontsize=20)
    curr_ax.set_ylabel('UMAP2',fontsize=20)
    curr_ax.legend(bbox_to_anchor=(1,.9),title='Orders',prop={'size': 12})#this is for the order umap
    plt.tight_layout()

#hyperparameter tuning https://umap-learn.readthedocs.io/en/latest/parameters.html
#function for performing grid search
def grid_search(n_neighbors,min_dist):
    X_scaled=StandardScaler().fit_transform(X)
    reducer=umap.UMAP(
            random_state=0,
            min_dist=min_dist,
            n_neighbors=n_neighbors
        )
    embedding=reducer.fit_transform(X_scaled,y=features['order_label'])
    features['x'],features['y']=embedding[:,0],embedding[:,1]
    fig,axs=plt.subplots(1,1,figsize=(8.5,6.5))
    plotter(axs)
    plt.savefig(f'../results/{n_neighbors}_nn_{min_dist}_md.png',dpi=400,bbox_inches='tight')

#parameters
#n_neighbors,min_dists=[10,50,200],[0.1,0.5,0.9]

#for i in n_neighbors:
#    for j in min_dists:
#        grid_search(i,j)

grid_search(400,.9)

features[['order','x','y']].to_csv('../data/order_umap.csv')