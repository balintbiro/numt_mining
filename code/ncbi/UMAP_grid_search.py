#import dependencies
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

#read in dataframe
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#add sizes
numts['gsize_comp']=numts['genomic_length']-(numts['gDNA_size (Mb)']*1000000)
numts['mtsize_comp']=numts['mitochondrial_length']-(numts['mtDNA_size (Mb)']*1000000)

#get taxonomy dictionaries
order_dict=pd.read_csv('../data/order_dict.csv',index_col=0)
genus_dict=pd.read_csv('../data/genus_dict.csv',index_col=0)
family_dict=pd.read_csv('../data/family_dict.csv',index_col=0)

#create filters to filter out missing taxonomy labels (it may have cause issues)
order_fil=numts['order_label']!=order_dict.loc[np.nan].values[0]
genus_fil=numts['genus_label']!=genus_dict.loc[np.nan].values[0]
family_fil=numts['family_label']!=family_dict.loc[np.nan].values[0]

numts=numts[order_fil][family_fil][genus_fil]

#create datasets
X_labeled=numts[[
    'score','eg2_value','e_value',#alignment scores
    'genomic_start','genomic_length','mitochondrial_length','genomic_size',#sequences features
    'numt_GC','upstream_GC','downstream_GC',#GCs
    'modk2','transversions','transitions',#pairwise divergence
    'uSW_mean', 'uSW_median', 'uRMs_count', 'uRMs_lengths',#upstream flanking features
    'dSW_mean', 'dSW_median', 'dRMs_count', 'dRMs_lengths',#downstream_flanking features
    'gnumt_relGC', 'u_relGC', 'd_relGC', 'grel_numt_size', 'mtrel_numt_size', 'mtnumt_relGC',#genomic data
    'u_1st_repeatl', 'u_2nd_repeatl', 'u_3rd_repeatl',
    'u_4th_repeatl', 'u_5th_repeatl','u_1st_repeatclassl','u_2nd_repeatclassl',#RM frequencies
    'gsize_comp','mtsize_comp',
    'gDNA_size (Mb)','mtDNA_size (Mb)',
        'genus_label','family_label','order_label','label']].dropna()

X=X_labeled.loc[:,X_labeled.columns!='order_label']

#define function for plotting the result
def plotter(coloring_label,curr_ax,title=None):
    kwargs={'edgecolor':'face'}
    scatter_plot=curr_ax.scatter(
        x='x',
        y='y',
        c=coloring_label,
        data=X,
        #cmap='Paired',
        s=5,
        **kwargs
    )
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
    embedding=reducer.fit_transform(X_scaled,y=X_labeled['order_label'])
    X['x']=embedding[:,0]
    X['y']=embedding[:,1]
    X['order_label']=X_labeled['order_label']
    X['genus_label']=X_labeled['genus_label']
    X['family_label']=X_labeled['family_label']
    X['label']=X_labeled['label']
    output=X[['order_label','genus_label','family_label','label','x','y']]
    output.to_csv('../data/coordinates.csv')
    fig,axs=plt.subplots(1,1,figsize=(10,10))
    plotter('order_label',axs)
    plt.savefig(f'../results/{n_neighbors}_nn_{min_dist}_md.png',dpi=400,bbox_inches='tight')

grid_search(400,0.99)
