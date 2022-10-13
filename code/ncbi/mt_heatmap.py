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

#get numt coverage
def get_numt_coverage(organism_name):
	numt_range=numt_ranges[organism_name]#all the nucleotides that are responsible for numtogenesis (redundant!)
	numt_coverage=[]
	#handling gaps and nucleotides
	nuc_counter=0
	for index in np.arange(0,len(msa.loc[organism_name])):
		base=msa.loc[organism_name][index]
		if base!='-':
			numt_coverage.append(numt_range.count(nuc_counter))
			nuc_counter+=1
		else:
			numt_coverage.append(0)
	return numt_coverage

#calculate numt coverages
numt_coverages=organism_names.apply(get_numt_coverage)
numt_coverages.index=organism_names
numt_coverages.to_csv('../results/mt_heatmap.csv')