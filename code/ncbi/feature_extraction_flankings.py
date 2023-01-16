#import dependencies
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import iFeatureOmegaCLI
from subprocess import call

#read in numts and random sequences
numts=pd.read_csv('../data/ncbi_numts_p26.csv')
randoms=pd.read_csv('../data/ml_input_randoms.csv')

#create filters; filter the lines that have at least 200 bp flanking regions
flanking=200
numt_upfil,numt_downfil,random_upfil,random_downfil=numts['upstream_5kb'].apply(lambda seq: len(seq)>flanking if type(seq)!= float else False),numts['downstream_5kb'].apply(lambda seq: len(seq)>flanking),randoms['upstream_size']>flanking,randoms['downstream_size']>flanking

#do the actual filtering; and sample random sequences so we have the same number of data points in both sets
numts=numts[numt_upfil][numt_downfil]
randoms=randoms[random_upfil][random_downfil].sample(n=len(numts),replace=False)

#add 200 bp flankings to each dataset (random and numt)
numts['u_flanking'],numts['d_flanking']=numts['upstream_5kb'].apply(lambda seq: seq.upper()[-200:]),numts['downstream_5kb'].apply(lambda seq: seq.upper()[:200])
randoms['u_flanking'],randoms['d_flanking']=randoms.apply(lambda row: row['sequence'][(row['upstream_size']-200):row['upstream_size']].upper(),axis=1),randoms.apply(lambda row: row['sequence'][(row['upstream_size']+row['sample_size']):(row['upstream_size']+row['sample_size']+200)].upper(),axis=1)

#add labels
numts['label'],randoms['label']=len(numts)*[1],len(randoms)*[0]

#reduce column numbers
numts,randoms=numts[['genomic_id','u_flanking','d_flanking','label']],randoms[['genomic_id','u_flanking','d_flanking','label']]

#merge the two dataframes
sequences=pd.concat([numts,randoms])

#function to write the sequences for further analysis
def write_seq(row,outfile):
	g_id=row['genomic_id']
	u_seq,d_seq=row['u_flanking'],row['d_flanking']
	if row['label']==1:
		outfile.write(f'>nu_{g_id}\n{u_seq}\n')
		outfile.write(f'>nd_{g_id}\n{d_seq}\n')
	else:
		outfile.write(f'>ru_{g_id}\n{u_seq}\n')
		outfile.write(f'>rd_{g_id}\n{d_seq}\n')

#create file for the sequences
with open('../data/flanking_sequences.fasta','w')as outfile:
	sequences.apply(write_seq,args=(outfile,),axis=1)

#create folder for descriptor dataframes
if os.path.exists('../data/flanking_features/')==False:
	os.mkdir('../data/flanking_features/')

#define descriptor types
descriptors=pd.Series([
	'Mismatch'
	])

#function for getting descriptors
def get_desc(descriptor):
	try:
		#create DNA instance with iFeatureOmegaCLI for numts
		dna=iFeatureOmegaCLI.iDNA('../data/flanking_sequences.fasta')
		#get parameters
		dna.import_parameters('../data/DNA_parameters_setting.json')
		dna.get_descriptor(descriptor)
		if dna.encodings is not None:
			dna.to_csv(f"../data/flanking_features/{descriptor.replace(' ','')}.csv",index=False,header=True)
	except:
		print(f'Could not generate {descriptor} related features!')

#generate descriptors
descriptors.apply(get_desc)

#get fasta headers aka genomic ids
records=list(SeqIO.parse('../data/flanking_sequences.fasta','fasta'))
genomic_ids=pd.Series(records).apply(lambda record: record.description)

files=pd.Series(['NAC.csv','Mismatch.csv','Kmertype1.csv','Kmertype2.csv','NMBroto.csv','Z_curve_9bit.csv','RCKmertype1.csv','RCKmertype2.csv',])

#get features
features=pd.concat(files.apply(lambda filename: pd.read_csv(f'../data/flanking_features/{filename}')).tolist(),axis=1)

#set fasta headers as indices
features.index=genomic_ids

#write features
features.to_csv('../data/flanking_features.csv')