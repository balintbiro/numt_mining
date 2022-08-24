#import dependencies
import numpy as np
import pandas as pd

#read in numts
numts=pd.read_csv('../data/ncbi_numts_p15.csv')

#read nuclear information
gDNA=pd.read_csv('../data/eukaryotes.txt',sep='\t',index_col=0)

#read mitochondrial information
mtDNA=mtdnas=pd.read_csv('../data/organelles.csv',index_col=0)

#function for getting scaffold counts
def gDNA_scaffolds(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		gDNA_feature=gDNA.loc[organism_name]['Scaffolds'].tolist()
		if type(gDNA_feature)==list:
			return np.mean(gDNA_feature)
		else:
			return gDNA_feature
	except:
		return np.nan

#add scaffold counts to numts df
numts['scaffolds']=numts['organism_name'].apply(gDNA_scaffolds)

#function for getting gDNA GC
def gDNA_GC(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		gDNA_feature=gDNA.loc[organism_name]['GC%']
		if type(gDNA_feature)==str:
			return float(gDNA_feature)
		else:
			gDNA_feature=gDNA.loc[organism_name]['GC%'].dropna()
			fil=gDNA_feature.apply(lambda item: item!='-')
			return np.mean(list(map(float,gDNA_feature[fil].tolist())))
	except:
		return np.nan

#add total gDNA GSc to numts df
numts['gDNA_GC']=numts['organism_name'].apply(gDNA_GC)

#function for getting gDNA sizes
def gDNA_size(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		gDNA_feature=gDNA.loc[organism_name]['Size (Mb)']
		if type(gDNA_feature)==float:
			return gDNA_feature
		else:
			return np.mean(gDNA_feature)
	except:
		return np.nan

#add gDNA_sizes to numts df
numts['gDNA_size (Mb)']=numts['organism_name'].apply(gDNA_size)

#function for getting mtDNA GC
def mtDNA_GC(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		mtDNA_feature=mtDNA.loc[organism_name]['GC%']
		if type(mtDNA_feature)==pd.Series:
			return np.mean(mtDNA_feature)
		else:
			return mtDNA_feature
	except:
		return np.nan

#add total gDNA GSc to numts df
numts['mtDNA_GC']=numts['organism_name'].apply(mtDNA_GC)

#function for getting mtDNA sizes
def mtDNA_size(organism_name):
	try:
		organism_name=organism_name[0].upper()+organism_name[1:].replace('_',' ')
		mtDNA_feature=mtDNA.loc[organism_name]['Size(Mb)']
		if type(mtDNA_feature)==float:
			return mtDNA_feature
		else:
			return np.mean(mtDNA_feature)
	except:
		return np.nan

#add mtDNA sizes to numts df
numts['mtDNA_size (Mb)']=numts['organism_name'].apply(mtDNA_size)

#add numts relative GC (compared to gDNA GC)
numts['gnumt_relGC']=numts['numt_GC']/(numts['gDNA_GC']/100)

#add upstream relative GC
numts['u_relGC']=numts['upstream_GC']/(numts['gDNA_GC']/100)

#add downstream rel GC
numts['d_relGC']=numts['downstream_GC']/(numts['gDNA_GC']/100)

#add rel numt size compared to genome size
numts['grel_numt_size']=numts['genomic_length']/(numts['gDNA_size (Mb)']*1000000)

#add rel numt size compared to mt size
numts['mtrel_numt_size']=numts['genomic_length']/(numts['mtDNA_size (Mb)']*1000000)

#add numts relative GC comapred to mt
numts['mtnumt_relGC']=numts['numt_GC']/(numts['mtDNA_GC']/100)

#wwrite output
numts.to_csv('../data/ncbi_numts_p16.csv',index=False)