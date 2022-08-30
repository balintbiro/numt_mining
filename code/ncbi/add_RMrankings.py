#import dependencies
import numpy as np
import pandas as pd

#reading in merged numts df
numts=pd.read_csv('../data/ncbi_numts_p16.csv')

#reading in the repeat files
upstream_repeats=pd.read_csv('../results/upstream_RM.csv')
downstream_repeats=pd.read_csv('../results/downstream_RM.csv')

#add flanking types
upstream_repeats['flanking']=len(upstream_repeats)*['upstream']
downstream_repeats['flanking']=len(downstream_repeats)*['downstream']

#merge the repeat dfs into one df
repeats=pd.concat([upstream_repeats,downstream_repeats])

#function for adding different levels of repeats
def add_RMranking(row,sequence_type):
	try:
		query_name=f"{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(len(row[sequence_type]))}_{str(int(row['mitochondrial_start']))}"
		subdf=repeats.loc[repeats['query_name']==query_name]
		if len(subdf)>0:
			urr=list(subdf['repeat'].value_counts().index[:5])#urr stands for upstream repeat ranks
			urcr=list(subdf['repeat_class'].value_counts().index[:5])#urcr stands for upstream repeat class ranks
			if len(urr)!=5:
				urr+=(5-len(urr))*[np.nan]
			elif len(urcr)!=5:
				urcr+=(5-len(urcr))*[np.nan]
			repeatranks=urr+urcr
			return repeatranks
		else:
			return 10*[np.nan]
	except:
		return 10*[np.nan]

#get upstream repeat and repeat class ranks
urrs=pd.DataFrame(numts.apply(add_RMranking,args=('upstream_5kb',),axis=1).tolist())#urs stands for upstream repeat ranks
urrs.columns=[
		'u_1st_repeat','u_2nd_repeat','u_3rd_repeat','u_4th_repeat','u_5th_repeat',
'u_1st_repeatclass','u_2nd_repeatclass','u_3rd_repeatclass','u_4th_repeatclass','u_5th_repeatclass'
]

#get downstream repeat and repeat class ranks
drrs=pd.DataFrame(numts.apply(add_RMranking,args=('downstream_5kb',),axis=1).tolist())#urs stands for upstream repeat ranks
drrs.columns=[
		'd_1st_repeat','d_2nd_repeat','d_3rd_repeat','d_4th_repeat','d_5th_repeat',
'd_1st_repeatclass','d_2nd_repeatclass','d_3rd_repeatclass','d_4th_repeatclass','d_5th_repeatclass'
]

#get merged repeatranks
repeatranks=pd.concat([urrs,drrs],axis=1)

#function for creating labels for each repeat column
def get_labels(column_name):
	rm_dict=pd.Series(
			data=np.arange(0,len(repeatranks[column_name].unique())),
			index=repeatranks[column_name].unique()
		)
	return list(rm_dict[repeatranks[column_name]].values)

#create RM rank labels
repeatrank_labels=pd.DataFrame(pd.Series(repeatranks.columns).apply(get_labels).tolist()).T
repeatrank_labels.columns=[
	'u_1st_repeatl','u_2nd_repeatl','u_3rd_repeatl','u_4th_repeatl','u_5th_repeatl',
'u_1st_repeatclassl','u_2nd_repeatclassl','u_3rd_repeatclassl','u_4th_repeatclassl','u_5th_repeatclassl',
'd_1st_repeatl','d_2nd_repeatl','d_3rd_repeatl','d_4th_repeatl','d_5th_repeatl',
'd_1st_repeatclassl','d_2nd_repeatclassl','d_3rd_repeatclassl','d_4th_repeatclassl','d_5th_repeatclassl'
]

#merge together numts and repeatrank labels
numts=pd.concat([numts,repeatrank_labels],axis=1)

#write_output
numts.to_csv('../data/ncbi_numts_p26.csv',index=False)
repeatranks.to_csv('../results/repeatranks.csv',index=False)