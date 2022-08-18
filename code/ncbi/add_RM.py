#import dependencies
import os
import pandas as pd
from subprocess import call

#reading in merged numts df
numts=pd.read_csv('../data/ncbi_numts.csv')

#create folder for RepeatMasker files
if os.path.isdir('../data/RM_files/')==False:
    os.mkdir('../data/RM_files/')

#function for preparing the output for RepeatMasker analysis, run RepeatMasker itself and process outputs
def repeatmasker(organism_name,sequence_type):
	global RM_dfs
	try:
		subdf=numts.loc[numts['organism_name']==organism_name]
		#write the sequences into a fasta file
		with open('../data/RM_files/RM_input.fa','w')as outfile:
			subdf.apply(
					lambda row: outfile.write(
							f">{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(len(row[sequence_type]))}_{str(int(row['mitochondrial_start']))}\n{row[sequence_type]}\n" #header consists of genomic id+flanking length, genomic and mitochondrial starts (to make sure of uniqueness)
						)
					if type(row[sequence_type])==str #in some cases there is no flanking region and so np.nans are present
					else print(f"{organism_name} doesnt have {sequence_type} under the genomic ID of {row['genomic_id']} in the position of {str(row['genomic_start']).upper()}"),
					axis=1
				)
		#change directory and call RepeatMasker
		call(f"""cd ../data/RM_files; RepeatMasker -species {organism_name.replace('_',' ')} RM_input.fa""",shell=True)
		#transform RepeatMasker.out file into a tab separated file for better handling
		call("""cat ../data/RM_files/RM_input.fa.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' > ../data/RM_files/RM_input.tab""",shell=True)
		#add columns https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
		repeats=pd.read_csv('../data/RM_files/RM_input.tab',sep='\t',skiprows=2,header=None)
		#add organism name to repeats df
		repeats['organism_name']=len(repeats)*[organism_name]
		#append individual repeat df to a common df list
		RM_dfs.append(repeats)
		#clear the directory
		call('rm -r ../data/RM_files/*',shell=True)
	except:
		pass

#get organism names
organism_names=pd.Series(os.listdir('../data/alignments/')).apply(lambda name: '_'.join(name.split('_')[:2]).lower())

######################################################################################################
#create global variable for upstream RM dfs
RM_dfs=[]
#make the function work for upstream repeats
organism_names.apply(repeatmasker,args=('upstream_5kb',))
#create a merged df for upstream repeats
upstream_repeats=pd.concat(RM_dfs)
#add column names for upstream repeats
upstream_repeats.columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
#write upstream outputs
upstream_repeats.to_csv('../results/upstream_RM.csv',index=False)

######################################################################################################
#create global variable for downstream RM dfs
RM_dfs=[]
#make the function work for downstream repeats
organism_names.apply(repeatmasker,args=('downstream_5kb',))
#create a merged df for downstream repeats
downstream_repeats=pd.concat(RM_dfs)
#add column names for upstream repeats
downstream_repeats.columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
#write upstream outputs
downstream_repeats.to_csv('../results/downstream_RM.csv',index=False)