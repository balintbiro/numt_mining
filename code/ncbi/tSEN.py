import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

def read_csv(filename):
    global indicer
    try:
        df=pd.read_csv(f'../data/alignments/{filename}')
        df['label']=len(df)*[indicer]
        indicer+=1
        df['GC']=df['genomic_sequence'].apply(lambda seq:(seq.upper().count('G')+seq.upper().count('C'))/len(seq))

        return df
    except:
        return np.nan

indicer=0
dfs=pd.Series(os.listdir('../data/alignments/')).apply(read_csv)

dfs=dfs.dropna()
dfs

df=pd.concat(dfs.tolist())

X=df[['score','eg2_value','e_value','genomic_start','genomic_length','genomic_size','GC']]
y=df['label']

plt.style.use('fivethirtyeight')

def tSNE(indexer):
	X_normalized=StandardScaler().fit(X).transform(X)
	tsne = TSNE(
		random_state = 0,
		perplexity=perplexity_values[indexer],
		learning_rate=learning_rate_values[indexer],
		n_iter=n_iter_values[indexer],
		n_jobs=-1
		)
	X_tsne = tsne.fit_transform(X_normalized)
	X['x']=X_tsne[:,0]
	X['y']=X_tsne[:,1]
	X['label']=y
	fig,axs=plt.subplots(1,1,figsize=(8,8))
	sns.scatterplot(
	    x='x',
	    y='y',
	    hue='label',
	    data=X,
	    palette='Dark2',
	    alpha=.7,
	    ax=axs
	)
	plt.legend([])
	plt.savefig(f'../results/tsnes/pp{perplexity_values[indexer]}_ll{learning_rate_values[indexer]}_ni{n_iter_values[indexer]}.png',dpi=300)

perplexity_values=np.linspace(5,50,dtype=int)
learning_rate_values=np.linspace(10,1000,num=50,dtype=int)
n_iter_values=np.linspace(250,10000,num=50,dtype=int)

pd.Series(np.linspace(0,49,dtype=int)).apply(tSNE)