#import dependencies
import os
import numpy as np
import pandas as pd

#read in numts
numts=pd.read_csv('../data/ncbi_numts_p14.csv')

#get labels based on organism name
indices=pd.Series(
	data=np.arange(0,len(numts['organism_name'].unique())),
	index=numts['organism_name'].unique()
	)

#add labels to numts df
numts['label']=numts['organism_name'].apply(lambda organism_name: indices[organism_name])

#read in phylogenetic ranks
rankings=pd.read_csv('../data/rankedlineage.dmp',sep='\t',header=None)
fil=rankings.iloc[0].apply(lambda item:item!='|')#filter ranking columns
rankings=rankings[rankings.columns[fil]]#filter ranking columns
rankings.columns=['tax_id','Tax_name','species','genus','family','order','class','phylum','kingdom','superkingdom']#add ranking columns
rankings=rankings.set_index('tax_id')#set tax ids as indices

#read in Taxids
tax_ids=pd.read_csv('../data/eukaryotes.txt',sep='\t',index_col=0)

#get labels based on organism name
indices=pd.Series(
	data=np.arange(0,len(numts['organism_name'].unique())),
	index=numts['organism_name'].unique()
	)

#add labels to numts df
numts['label']=numts['organism_name'].apply(lambda organism_name: indices[organism_name])

#function for getting tax ids for every organism names in numts df
def get_taxid(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		tax_id=tax_ids.loc[organism_name]['TaxID'].tolist()
		if type(tax_id)==int:
			return tax_id
		else:
			return tax_id[0]
	except:
		return np.nan

#get tax_ids
taxids=numts['organism_name'].apply(get_taxid)

#function for getting phylogenetic ranks for each tax id
def get_phylorank(tax_id):
	try:
		return rankings.loc[tax_id].tolist()
	except:
		return len(rankings.columns)*[np.nan]

#get phylogenetic ranks based on previously determined tax ids
phyloranks=taxids.apply(get_phylorank)

#create dataframe from phylogenetic ranks
phyloranks=pd.DataFrame(
		data=phyloranks.tolist(),
		columns=rankings.columns
	)

#merge the two dfs
numts=pd.concat([numts,phyloranks],axis=1)

#delete unnecessary columns
numts=numts.drop(['Unnamed: 0', 'Unnamed: 0.1','species','class','phylum','kingdom','superkingdom','Tax_name'],axis=1)

#function for creating dictionaries for each phylo column to code them to number
def create_dict(column_name):
	return pd.Series(
			data=np.arange(0,len(numts[column_name].unique())),
			index=numts[column_name].unique()
		)

#create phylo dictionaries
genus_dict=create_dict('genus')
family_dict=create_dict('family')
order_dict=create_dict('order')

#translate each phylo type to numbers and add them to numts df
numts['genus_label']=genus_dict[numts['genus']].values
numts['family_label']=family_dict[numts['family']].values
numts['order_label']=order_dict[numts['order']].values

#write output
numts.to_csv('../data/ncbi_numts_p15.csv',index=False)