#import dependencies
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, LabelEncoder

#read in dataframe
features=pd.read_csv('../data/flanking_features.csv',index_col=0)

#add label names
def get_label(name):
    if name[0]=='n':
        if name[1]=='u':
            return 'NUMT upstream'
        else:
            return 'NUMT downstream'
    else:
        if name[1]=='u':
            return 'random upstream'
        else:
            return 'random downstream'
    
features['label_name']=pd.Series(features.index).apply(get_label).values

#create labels and add them to the dataframe
features['label']=LabelEncoder().fit_transform(features['label_name'].values)

X=features.drop(['label','label_name'],axis=1)

#define function for plotting the result
def plotter(curr_ax,title=None):
    kwargs={'edgecolor':'face'}
    scplot=sns.scatterplot(x='x',y='y',hue='label',data=features,s=10,palette='tab20',alpha=0.25,**kwargs)
    curr_ax.set(xticklabels=[],yticklabels=[])
    #scplot.legend(title='Sequence',fontsize=15,title_fontsize=15)#this is for the numt vs random umap
    curr_ax.set_xlabel('UMAP1',fontsize=20)
    curr_ax.set_ylabel('UMAP2',fontsize=20)
    curr_ax.legend(bbox_to_anchor=(1,.9),title='Flanking sequence',prop={'size': 12})#this is for the order umap
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
    embedding=reducer.fit_transform(X_scaled,y=features['label'])
    features['x'],features['y']=embedding[:,0],embedding[:,1]
    fig,axs=plt.subplots(1,1,figsize=(8.5,6.5))
    plotter(axs)
    plt.savefig(f'../results/flanking_{n_neighbors}_nn_{min_dist}_md.png',dpi=400,bbox_inches='tight')

#parameters
#n_neighbors,min_dists=[10,50,200],[0.1,0.5,0.9]

#for i in n_neighbors:
#    for j in min_dists:
#        grid_search(i,j)

grid_search(400,.9)

features[['label_name','x','y']].to_csv('../data/flanking_umap.csv')

####################################################################################################
#                                       visualisation                                              #
#           for just the visualization, copy the code from here plus the dependencies              #
####################################################################################################

features=pd.read_csv('../data/flanking_umap.csv',index_col=0)
features['2label_name']=features['label_name'].apply(lambda name: name.rsplit()[0])

fig,axs=plt.subplots(1,2,figsize=(12.5,6.),sharey=True)

kwargs={'edgecolor':'face'}
scplot1=sns.scatterplot(x='x',y='y',hue='label_name',data=features,s=10,ax=axs[0],palette='tab10',alpha=0.25,**kwargs)
axs[0].set(xticklabels=[],yticklabels=[])
axs[0].set_xlabel('UMAP1',fontsize=20)
axs[0].set_ylabel('UMAP2',fontsize=20)
scplot1.legend(title='Flanking sequence',prop={'size': 12})

scplot2=sns.scatterplot(x='x',y='y',hue='2label_name',data=features,s=10,ax=axs[1],palette='tab10',alpha=0.25,**kwargs)
axs[1].set_xlabel('UMAP1',fontsize=20)
axs[1].set(xticklabels=[],yticklabels=[])
scplot2.legend(title='Flanking sequence',prop={'size': 12})
plt.tight_layout()
plt.savefig('../results/flanking_400_nn_09_md.png',dpi=400)