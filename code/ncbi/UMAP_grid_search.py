#import dependencies
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

#read in dataframe
numts=pd.read_csv('../data/ncbi_numts_p14.csv',index_col=0)

#get labels based on organism name
indices=pd.Series(
	data=np.arange(0,len(numts['organism_name'].unique())),
	index=numts['organism_name'].unique()
	)

#add labels to numts df
numts['label']=numts['organism_name'].apply(lambda organism_name: indices[organism_name])

#dropnas since its not suitable
numts=numts.dropna()

#create datasets
X=numts[[
    'score','eg2_value','e_value',#alignment scores
    'genomic_start','genomic_length','mitochondrial_length','genomic_size',#sequences features
    'numt_GC','upstream_GC','downstream_GC',#GCs
    'modk2','transversions','transitions',#pairwise divergence
    'uSW_mean', 'uSW_median', 'uRMs_count', 'uRMs_lengths',#upstream flanking features
    'dSW_mean', 'dSW_median', 'dRMs_count', 'dRMs_lengths'#downstream_flanking features
        ]]
y=numts['label']

#conditional creation of UMAP results folder
if os.path.exists('../results/UMAPs/')==False:
    os.mkdir('../results/UMAPs/')

#function for performing grid search
def grid_search(neighbor_value,component_value):
    X_scaled=StandardScaler().fit(X).transform(X)
    reducer=umap.UMAP(
            random_state=0,
            n_neighbors=neighbor_value,
            n_components=component_value
        )
    embedding=reducer.fit_transform(X_scaled)
    plt.style.use('fivethirtyeight')
    fig,axs=plt.subplots(1,1,figsize=(10,10))
    axs.scatter(
            embedding[:, 0],
            embedding[:, 1],
            c=y,
            cmap='magma'
        )
    axs.axis('off')
    plt.legend().remove()
    plt.tight_layout()
    plt.savefig(f'../results/UMAPs/{neighbor_value}nn_{component_value}nc.png',dpi=450)

#apply grid_search function for different paramters
#for neighbor_value in np.linspace(2,100,5,dtype=int):
#    for component_value in np.linspace(2,100,5,dtype=int):
#        grid_search(neighbor_value,component_value)
grid_search(2,2)