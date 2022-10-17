#import dependencies
import numpy as np
import pandas as pd
from Bio import AlignIO

#read in numts
numts=pd.read_csv('../data/ncbi_numts_p26.csv')
numts=numts[['organism_name','mitochondrial_start','mitochondrial_length']]

#calculate numt ranges
numts['numt_ranges']=numts.apply(
		lambda row: list(np.arange(row['mitochondrial_start'],row['mitochondrial_start']+row['mitochondrial_length'])),
		axis=1
	)

#get organism names
organism_names=pd.Series(numts['organism_name'].unique())

#function for concatenating numt ranges for each organism
def get_numt_ranges(organism_name):
	subdf=numts.loc[numts['organism_name']==organism_name]
	return sum(subdf['numt_ranges'].tolist(),[])

#make the function work
numt_ranges=organism_names.apply(get_numt_ranges)
numt_ranges.index=organism_names

#read in alignments
msa=AlignIO.read(
		open('../results/aligned_mtDNAs.fa'),
		format='fasta'
	)

#add ids to msa
msa_ids=[]
for record in msa:
	msa_ids.append(record.id.lower())
msa=pd.DataFrame(msa)
msa.index=msa_ids

#transform msa-do not show individual nucleotides but their position in mitochondrion
def transform_msa(organism_name):
	alignment_range=msa.loc[organism_name]
	alignment_range[alignment_range!='-']=np.arange(0,len(alignment_range[alignment_range!='-']))
	msa.loc[organism_name]=alignment_range

#function for getting numt_coverages
def get_numt_coverage(organism_name,output_file):
	numt_range=pd.Series(numt_ranges[organism_name]).value_counts()
	numt_range['-']=-1#add minus 1 for gaps
	msa_range=msa.loc[organism_name]
	numt_coverage=list(numt_range.reindex(msa_range.values).fillna(0).values)
	numt_coverage=[organism_name]+numt_coverage
	numt_coverage=pd.DataFrame([numt_coverage])
	numt_coverage.to_csv(output_file,mode='a',index=False,header=False)

#transform msa from alignment to mt positions
organism_names.apply(transform_msa)

#create empty df
heatmaps=pd.DataFrame(data=[],columns=np.arange(0,msa.shape[1]))
#write it to csv
heatmaps.to_csv('../results/heatmap.csv',index=False)

#transform df to numt coverages
organism_names.apply(get_numt_coverage,args=('../results/heatmap.csv',))
