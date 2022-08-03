#import dependencies
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

#read in dataframe
#import numts
numts=pd.read_csv('../data/ncbi_numts_modk2_added.csv',index_col=0)
numts['score']=numts.index
numts.index=np.arange(0,len(numts))

#create datasets
X=numts[[
    'score','eg2_value','e_value','genomic_start','genomic_length','mitochondrial_length','GC'
        ]]#
y=numts['label']

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
	axs.axis('off')
	plt.legend().remove()
	plt.tight_layout()
	plt.savefig(f'../results/tsnes/pp{perplexity_value}_without_modk2.png',dpi=450)

#apply function
pd.Series(np.linspace(5,100,10,dtype=int)).apply(grid_search)