#import dependencies
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

#read in dataframe
numts=pd.read_csv('../results/tsne_input.csv',index_col=0)

X=numts[['score','eg2_value','e_value','genomic_start','genomic_length','mitochondrial_length','numt_gc',
        'g_gc','full_g_size','scaffolds_number','genes_number','proteins_number','full_mt_size','mt_gc',
        'mt_cds']]
y=numts['organism_label']

#function for performing grid search
def grid_search(perplexity_value):
	X_normalized=StandardScaler().fit(X).transform(X)
	tsne = TSNE(random_state = 0,perplexity=perplexity_value,n_jobs=-1)
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
	plt.legend([])
	plt.tight_layout()
	plt.savefig(f'../results/tsnes/pp{perplexity_value}.png',dpi=450)

#apply function
pd.Series(np.arange(5,50)).apply(grid_search)