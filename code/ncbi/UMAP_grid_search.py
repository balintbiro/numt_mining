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
    'u_1st_repeatl', 'u_2nd_repeatl', 'u_3rd_repeatl','u_4th_repeatl', 'u_5th_repeatl','u_1st_repeatclassl','u_2nd_repeatclassl',#RM frequencies
    'gsize_comp','mtsize_comp',
    'gDNA_size (Mb)','mtDNA_size (Mb)',
        'genus_label','family_label','order_label','label']]

#sample df
X_labeled=X_labeled.sample(frac=0.25)

X_labeled=X_labeled[[
    'score','eg2_value','e_value',#alignment scores
    'genomic_start','genomic_length','mitochondrial_length','genomic_size',#sequences features
    'numt_GC','upstream_GC','downstream_GC',#GCs
    'modk2','transversions','transitions',#pairwise divergence
    'uSW_mean', 'uSW_median', 'uRMs_count', 'uRMs_lengths',#upstream flanking features
    'dSW_mean', 'dSW_median', 'dRMs_count', 'dRMs_lengths',#downstream_flanking features
    'gnumt_relGC', 'u_relGC', 'd_relGC', 'grel_numt_size', 'mtrel_numt_size', 'mtnumt_relGC',#genomic data
    'u_1st_repeatl', 'u_2nd_repeatl', 'u_3rd_repeatl','u_1st_repeatclassl',#'u_4th_repeatl', 'u_5th_repeatl',#'u_2nd_repeatclassl',#RM frequencies
    'genus_label','family_label','order_label','label',
    'gsize_comp','mtsize_comp',
    'gDNA_size (Mb)','mtDNA_size (Mb)'
    ]].dropna()

X=X_labeled[[
    'score','eg2_value','e_value',#alignment scores
    'genomic_start','genomic_length','mitochondrial_length','genomic_size',#sequences features
    'numt_GC','upstream_GC','downstream_GC',#GCs
    'modk2','transversions','transitions',#pairwise divergence
    'uSW_mean', 'uSW_median', 'uRMs_count', 'uRMs_lengths',#upstream flanking features
    'dSW_mean', 'dSW_median', 'dRMs_count', 'dRMs_lengths',#downstream_flanking features
    'gnumt_relGC', 'u_relGC', 'd_relGC', 'grel_numt_size', 'mtrel_numt_size', 'mtnumt_relGC',#genomic data
    'u_1st_repeatl', 'u_2nd_repeatl', 'u_3rd_repeatl','u_1st_repeatclassl',#'u_4th_repeatl', 'u_5th_repeatl',#'u_2nd_repeatclassl',#RM frequencies
    'genus_label','label',#'family_label','order_label',
    'gsize_comp','mtsize_comp',
    'gDNA_size (Mb)','mtDNA_size (Mb)'
    ]]

#conditional creation of UMAP results folder
if os.path.exists('../results/UMAPs/')==False:
    os.mkdir('../results/UMAPs/')

#define function for plotting the result
def plotter(coloring_label,color_palette,curr_ax,title):
    sns.scatterplot(
        x='x',
        y='y',
        hue=coloring_label,
        data=X,
        palette=color_palette,
        alpha=.7,
        ax=curr_ax
    )
    curr_ax.set_title(title)
    curr_ax.axis('off')
    curr_ax.get_legend().remove()
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
    embedding=reducer.fit_transform(X_scaled)
    X['x']=embedding[:,0]
    X['y']=embedding[:,1]
    X['family_label']=X_labeled['family_label']
    X['order_label']=X_labeled['order_label']
    X['genus_label']=X_labeled['genus_label']
    X['species_label']=X_labeled['label']
    plt.style.use('fivethirtyeight')
    fig,axs=plt.subplots(2,2,figsize=(10,10))
    plotter('family_label','tab20',axs[0,0],'family')
    plotter('order_label','tab20_r',axs[0,1],'order')
    plotter('genus_label','Paired',axs[1,0],'genus')
    plotter('species_label','Paired_r',axs[1,1],'species')
    plt.savefig(f'../results/UMAPs/{n_neighbors}_nn_{min_dist}_md.png',dpi=450)

#for n_neighbor in [50,100,200]:
#    for min_dist in [0.8,0.99]:
#        grid_search(n_neighbor,min_dist)

grid_search(int(len(X)/4),0.8)