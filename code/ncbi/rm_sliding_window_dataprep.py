#import dependencies
import os
import numpy as np
import pandas as pd
from subprocess import call

#reading in merged numts df
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#create function for getting longest chr
def longest_chr(organism_name):
	fil=numts['organism_name']==organism_name
	longest_chr_id=numts[fil].sort_values(by='genomic_size',ascending=False)['genomic_id'].tolist()[0]
	return longest_chr_id

#desired organism names
organism_names=pd.Series([
		'sus_scrofa',
		'bos_taurus',
		'mus_musculus',
		'felis_catus',
		'rattus_norvegicus',
		'canis_lupus'
	])

#make the function work
longest_chr_ids=pd.Series(
	data=organism_names.apply(longest_chr)
	)
longest_chr_ids.index=organism_names

#create folder for RepeatMasker files
if os.path.isdir('../data/RM_files/')==False:
    os.mkdir('../data/RM_files/')

#function for preparing the output for RepeatMasker analysis, run RepeatMasker itself and process outputs
def repeatmasker(organism_name,output_file):
	longest_chr_id=longest_chr_ids[organism_name]
	subdf=numts.loc[numts['organism_name']==organism_name]
	subdf=subdf.loc[subdf['genomic_id']==longest_chr_id]
	#write the sequences into a fasta file
	#header consists of genomic id+flanking length, genomic and mitochondrial starts (to make sure of uniqueness)
	with open('../data/RM_files/RM_input.fa','w')as outfile:
		subdf.apply(
				lambda row: outfile.write(
						f">{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(int(row['mitochondrial_start']))}\n{row['upstream_5kb'].replace('-','')+row['genomic_sequence'].replace('-','')+row['downstream_5kb'].replace('-','')}\n"
					),
				axis=1
			)
	#change directory and call RepeatMasker
	call(f"""cd ../data/RM_files; RepeatMasker -species '{organism_name.replace('_',' ')}' RM_input.fa""",shell=True)
	#transform RepeatMasker.out file into a tab separated file for better handling
	call("""cat ../data/RM_files/RM_input.fa.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' > ../data/RM_files/RM_input.tab""",shell=True)
	#add columns https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
	repeats=pd.read_csv('../data/RM_files/RM_input.tab',sep='\t',skiprows=2,header=None)
	#add organism name to repeats df
	repeats['organism_name']=len(repeats)*[organism_name]
	#add current repeats to an existing output file
	repeats.to_csv(output_file,mode='a',index=False,header=False)
	#clear the directory
	call('rm -r ../data/RM_files/*',shell=True)

#create df for longest chr repeats
repeats=pd.DataFrame(
	data=[],
	columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
	)

#write repeats csv into a file; this df is going to be prolonged with the specific RM dfs
repeats.to_csv('../results/longest_chr_RM.csv',index=False)
#make the function work
organism_names.apply(repeatmasker,args=('../results/longest_chr_RM.csv',))