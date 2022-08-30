#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

#read in dataframe
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

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
        'genus_label','family_label','order_label','label']]

#dropnas since its not suitable for tsne
X_labeled=X_labeled.dropna()

X=X_labeled[[
    #'score','eg2_value','e_value',#alignment scores
    #'genomic_start','genomic_length','mitochondrial_length','genomic_size',#sequences features
    #'numt_GC','upstream_GC','downstream_GC',#GCs
    'modk2','transversions','transitions',#pairwise divergence
    'uSW_mean', 'uSW_median', 'uRMs_count', 'uRMs_lengths',#upstream flanking features
    'dSW_mean', 'dSW_median', 'dRMs_count', 'dRMs_lengths',#downstream_flanking features
    'gnumt_relGC', 'u_relGC', 'd_relGC', 'grel_numt_size', 'mtrel_numt_size', 'mtnumt_relGC',#genomic data
    'u_1st_repeatl', 'u_2nd_repeatl', 'u_3rd_repeatl','u_4th_repeatl', 'u_5th_repeatl','u_1st_repeatclassl','u_2nd_repeatclassl',#RM frequencies
    ]]

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
			init='pca',
			n_jobs=-1
		)
	X_tsne = tsne.fit_transform(X_normalized)
	X['x']=X_tsne[:,0]
	X['y']=X_tsne[:,1]
	X['family_label']=X_labeled['family_label']
	X['order_label']=X_labeled['order_label']
	X['genus_label']=X_labeled['genus_label']
	X['species_label']=X_labeled['label']
	plt.style.use('fivethirtyeight')
	fig,axs=plt.subplots(1,1,figsize=(10,10))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='family_label',
	    data=X,
	    palette='tab20',
	    alpha=.7,
	    ax=axs
	)
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tSNEs/family_{perplexity_value}pp_{learning_rate_value}lr.png',dpi=250)

	fig,axs=plt.subplots(1,1,figsize=(10,10))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='order_label',
	    data=X,
	    palette='tab20_r',
	    alpha=.7,
	    ax=axs
	)
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tSNEs/order_{perplexity_value}pp_{learning_rate_value}lr.png',dpi=250)

	fig,axs=plt.subplots(1,1,figsize=(10,10))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='genus_label',
	    data=X,
	    palette='Paired',
	    alpha=.7,
	    ax=axs
	)
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tSNEs/genus_{perplexity_value}pp_{learning_rate_value}lr.png',dpi=250)

	fig,axs=plt.subplots(1,1,figsize=(10,10))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='species_label',
	    data=X,
	    palette='Paired_r',
	    alpha=.7,
	    ax=axs
	)
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tSNEs/species_{perplexity_value}pp_{learning_rate_value}lr.png',dpi=250)

#apply function. Optimal hyperparameters: https://www.nature.com/articles/s41467-019-13056-x
for perplexity_value in np.linspace(5,len(numts)/100,4,dtype=int):
	for learning_rate_value in np.linspace(10,len(numts)/12,4,dtype=int):
		grid_search(perplexity_value,learning_rate_value)