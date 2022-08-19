#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
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

#dropnas since its not suitable for tsne
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

#conditional creation of tSNE results folder
if os.path.exists('../results/tSNEs/')==False:
    os.mkdir('../results/tSNEs/')

#function for performing grid search
def grid_search(perplexity_value,learning_rate_value):
	X_normalized=StandardScaler().fit(X).transform(X)
	tsne = TSNE(
			random_state = 0,
			perplexity=perplexity_value,
			learning_rate=learning_rate_value,
			n_jobs=-1
		)
	X_tsne = tsne.fit_transform(X_normalized)
	X['x']=X_tsne[:,0]
	X['y']=X_tsne[:,1]
	X['label']=y
	plt.style.use('fivethirtyeight')
	fig,axs=plt.subplots(1,1,figsize=(10,10))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='label',
	    data=X,
	    palette='Paired',
	    alpha=.7,
	    ax=axs
	)
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tSNEs/{perplexity_value}pp_{learning_rate_value}lr.png',dpi=450)

#apply function. Optimal hyperparameters: https://www.nature.com/articles/s41467-019-13056-x
for perplexity_value in np.linspace(5,len(numts)/100,5,dtype=int):
	for learning_rate_value in np.linspace(10,len(numts)/12,5,dtype=int):
		grid_search(perplexity_value,learning_rate_value)
#grid_search((len(numts)/100),(len(numts)/12))