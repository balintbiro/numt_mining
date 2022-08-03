#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

#read in csv
numts=pd.read_csv('../data/mice_numts.csv')

#create datasets
X=numts[[
    'score','eg2_value','e_value','g_start','g_length','mt_start','mt_length','GC','modk2','transitions','transversions'
        ]]#
y=numts['label']

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
	plt.savefig(f'../results/mice_tsnes/{perplexity_value}pp_{learning_rate_value}lr.png',dpi=450)

#apply function. Optimal hyperparameters: https://www.nature.com/articles/s41467-019-13056-x
for perplexity_value in np.linspace(5,len(numts)/100,7,dtype=int):
	for learning_rate_value in np.linspace(10,len(numts)/12,7,dtype=int):
		grid_search(perplexity_value,learning_rate_value)

grid_search(50,round(len(numts)/12,0))