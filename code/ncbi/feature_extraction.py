#import dependencies
import os
import pandas as pd
import iFeatureOmegaCLI
from subprocess import call

#read in numts and random sequences
numts=pd.read_csv('../data/ncbi_numts_p26.csv')
randoms=pd.read_csv('../data/ml_input_randoms.csv')

#get just the sequences and genomic ids
numt_sequences=numts[['genomic_id','genomic_sequence']]
random_sequences=randoms[['genomic_id','sequence']].sample(n=len(numt_sequences),replace=False)

#download the parameter file
if os.path.exists('../data/DNA_parameters_setting.json')==False:
	call('wget --directory-prefix=../data/ https://github.com/Superzchen/iFeatureOmega-CLI/tree/main/parameters/DNA_parameters_setting.json',shell=True)

#write just the sequences
with open('../data/numt_sequences.fasta','w')as outfile:
	numt_sequences.apply(lambda row: outfile.write(
			f">{row['genomic_id']}_{row.name}\n{row['genomic_sequence'].upper().replace('-','')}\n",
		),axis=1)

with open('../data/random_sequences.fasta','w')as outfile:
	random_sequences.apply(lambda row: outfile.write(
			f">{row['genomic_id']}_{row.name}\n{row['sequence'].upper().replace('-','')}\n",
		),axis=1)

#create folder for descriptor dataframes
if os.path.exists('../data/features/')==False:
	os.mkdir('../data/features/')

#define descriptor types already done 
#subsequence just eats way too much memory
descriptors=pd.Series([
	'Kmer type 1','Kmer type 2','RCKmer type 1','RCKmer type 2','Mismatch','NAC','Z_curve_9bit','Z_curve_12bit','Z_curve_36bit','Z_curve_48bit','Z_curve_144bit','CKSNAP type 1','CKSNAP type 2','MMI','NMBroto','ASDC',
		'ANF','ENAC','binary','PS2','PS3','PS4','NCP','EIIP','PseEIIP',
		'DBE','LPDF','DPCP','DPCP type2','TPCP','TPCP type2','Moran','Geary','DAC','DCC','DACC',
		'TAC','TCC','TACC','PseDNC','PseKNC','PCPseDNC','PCPseTNC','SCPseDNC','SCPseTNC','PSTNPss',
		'PSTNPds','KNN'
	])

#function for getting descriptors
def get_desc(descriptor,sequence_type):
	try:
		#create DNA instance with iFeatureOmegaCLI for numts
		if sequence_type=='numt':
			dna=iFeatureOmegaCLI.iDNA('../data/numt_sequences.fasta')
		else:
			dna=iFeatureOmegaCLI.iDNA('../data/random_sequences.fasta')
		#get parameters
		dna.import_parameters('../data/DNA_parameters_setting.json')
		dna.get_descriptor(descriptor)
		if dna.encodings is not None:
			dna.to_csv(f"../data/features/{sequence_type}_{descriptor.replace(' ','')}.csv",index=False,header=True)
	except:
		print(f'Could not generate {descriptor} related features!')

#make the function work
#descriptors.apply(get_desc,args=('numt',))
descriptors.apply(get_desc,args=('random',))