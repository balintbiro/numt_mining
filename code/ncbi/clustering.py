#import dependencies
import os
import umap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_mutual_info_score

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
        'genus_label','family_label','order','order_label','label']].dropna()

X_labeled=X_labeled.sample(n=int(len(X_labeled)/10))

#remove labels
X=X_labeled.loc[:,~X_labeled.columns.isin(['order','order_label'])]

#umap
numt_umap=umap.UMAP(random_state=0).fit_transform(X)

#normalize data
X_scaled=StandardScaler().fit_transform(X)

#pca
numt_pca=PCA(n_components=.95).fit_transform(X_scaled)

#get dbscan labels
dbscan_labels=DBSCAN(n_jobs=-1).fit_predict(numt_pca)

#get the datapoints that have been clustered
clustered=dbscan_labels>-1

kwargs={'edgecolor':'face'}

#visualize umap results
fig,axs=plt.subplots(1,1)
sns.scatterplot(x=numt_umap[:,0],y=numt_umap[:,1],hue=X_labeled['order'],ax=axs,palette='tab20',**kwargs)
plt.legend(bbox_to_anchor=(1.05,1),fontsize=10);plt.tight_layout()
fig.savefig('../results/UMAP.png',dpi=400)

#visualize PCA results
fig,axs=plt.subplots(1,1)
sns.scatterplot(x=numt_pca[:,0],y=numt_pca[:,1],hue=X_labeled['order'],ax=axs,palette='tab20',**kwargs)
plt.legend(bbox_to_anchor=(1.05,1),fontsize=10)
plt.tight_layout()
fig.savefig('../results/PCA.png',dpi=400)

#visualize DBSCAN clustered UMAP
fig,axs=plt.subplots(1,1)
sns.scatterplot(x=numt_umap[~clustered,0],y=numt_umap[~clustered,1],size=len(numt_umap[~clustered,1]),sizes=(10,20),hue=len(numt_umap[~clustered,1])*['black'],palette=['black'],legend=False,ax=axs,**kwargs)
sns.scatterplot(x=numt_umap[clustered,0],y=numt_umap[clustered,1],size=len(numt_umap[clustered,1]),sizes=(100,200),hue=X_labeled['order'][clustered],ax=axs,palette='tab20',**kwargs)
plt.legend(title='Order',bbox_to_anchor=(1.05,1),fontsize=10)
plt.tight_layout()
fig.savefig('../results/DBSCAN.png',dpi=400)

#grid search for UMAP and DBSCAN
min_distances=[0.1, 0.25, 0.5, 0.8, 0.99]
n_neighbors=np.linspace(555,1620,5,dtype=int)
epsilons=np.round(np.linspace(.25,20,10),2)
min_samples=np.linspace(2,15,10,dtype=int)

grid_search_results=pd.DataFrame(
		data=[],
		columns=['mn_dst','n_nghbrs','eps','mn_smpl','ami_cls','snr','ami']
	)

grid_search_results.to_csv('../results/DBSCAN_results.csv',index=False)

for min_distance in min_distances:
	for n_neighbor in n_neighbors:
		for eps in epsilons:
			for min_sample in min_samples:
				numt_umap=umap.UMAP(random_state=0,min_dist=min_distance,n_neighbors=n_neighbor).fit_transform(X)
				dbscan_labels=DBSCAN(eps=eps,min_samples=min_sample,n_jobs=-1).fit_predict(numt_pca)
				clustered=dbscan_labels>-1
				res1=adjusted_mutual_info_score(X_labeled['order'],dbscan_labels)
				res2=adjusted_mutual_info_score(X_labeled['order'][clustered],dbscan_labels[clustered])
				snr=len(X_labeled['order'][clustered])/len(X_labeled['order'])
				result=pd.DataFrame([[min_distance,n_neighbor,eps,min_sample,res2,snr,res1]])
				result.to_csv('../results/DBSCAN_results.csv',mode='a',index=False,header=False)


#import plotly.express as px
#fig=px.parallel_coordinates(
#		cv_res,
#		color='ami'
#	)
#fig.update_layout(
#		height=400,
#		font=dict(
#				family='Arial',
#				size=20,
#				color='#000000'
#			)
#	)

#fig.write_image(
#		'../results/DBSCAN_grid_search.png',
#		validate=True,
#		width=600,
#		height=400,
#		scale=3
#	)



