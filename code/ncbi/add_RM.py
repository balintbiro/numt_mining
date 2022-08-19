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
def repeatmasker(organism_name,sequence_type,output_file):
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
		#add current repeats to an existing output file
		repeats.to_csv(output_file,mode='a',index=False,header=False)
		#clear the directory
		call('rm -r ../data/RM_files/*',shell=True)
	except:
		pass

#get organism names
organism_names=pd.Series(os.listdir('../data/alignments/')).apply(lambda name: '_'.join(name.split('_')[:2]).lower())

######################################################################################################
#create a df for upstream repeats
#columns are from https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
upstream_repeats=pd.DataFrame(
	data=[],
	columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
	)
#write upstream csv into a file; this df is going to be prolonged with the species specific RM dfs
upstream_repeats.to_csv('../results/upstream_RM.csv',index=False)
#make the function work for upstream repeats
organism_names.apply(repeatmasker,args=('upstream_5kb','../results/upstream_RM.csv',))

######################################################################################################
#create a df for downstream repeats
downstream_repeats=pd.DataFrame(
	data=[],
	columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
	)
#write downstream csv into a file; this df is going to be prolonged with the species specific RM dfs
upstream_repeats.to_csv('../results/downstream_RM.csv',index=False)
#make the function work for downstream repeats
organism_names.apply(repeatmasker,args=('downstream_5kb','../results/downstream_RM.csv',))

#reading in the repeat files
upstream_repeats=pd.read_csv('../results/upstream_RM.csv')
downstream_repeats=pd.read_csv('../results/downstream_RM.csv')

#reading in numts
numts=pd.read_csv('../data/ncbi_numts_p6.csv')

#add flanking types
upstream_repeats['flanking']=len(upstream_repeats)*['upstream']
downstream_repeats['flanking']=len(downstream_repeats)*['downstream']

#merge the repeat dfs into one df
repeats=pd.concat([upstream_repeats,downstream_repeats])

#function for adding different features of RepeatMasker to the "master" df
def add_RM(row,sequence_type):
	try:
		query_name=f"{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(len(row[sequence_type]))}_{str(int(row['mitochondrial_start']))}"
		subdf=repeats.loc[repeats['query_name']==query_name]
		if len(subdf)!=0:
			SW_mean=np.mean(subdf['SW_score'])
			SW_median=np.median(subdf['SW_score'])
			RMs_count=len(subdf)
			RMs_lengths=sum(subdf['piq_end']-subdf['piq_begin'])
			return [SW_mean,SW_median,RMs_count,RMs_lengths]
		else:
			return [np.nan,np.nan,np.nan,np.nan]
	except:
		return [np.nan,np.nan,np.nan,np.nan]

#get RM features for upstream flanking
uRM_features=pd.DataFrame(numts.apply(add_RM,args=('upstream_5kb',),axis=1).tolist())
uRM_features.columns=['uSW_mean','uSW_median','uRMs_count','uRMs_lengths']

#get RM features for downstream flanking
dRM_features=pd.DataFrame(numts.apply(add_RM,args=('downstream_5kb',),axis=1).tolist())
dRM_features.columns=['dSW_mean','dSW_median','dRMs_count','dRMs_lengths']

#add RM features to numts df
numts=pd.concat([numts,uRM_features],axis=1)
numts=pd.concat([numts,dRM_features],axis=1)

#write output into a csv p14 means the number of features
numts.to_csv('../data/ncbi_numts_p14.csv',index=False)