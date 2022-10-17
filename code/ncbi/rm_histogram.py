#import dependencies
import os
import pandas as pd
from subprocess import call

#reading in merged numts df
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#filter df
numts=numts.dropna(subset=['upstream_5kb','downstream_5kb'])
numts=numts[['organism_name','genomic_id','genomic_start','upstream_5kb','genomic_length','downstream_5kb','mitochondrial_start','genomic_sequence']]
numts=numts[numts['upstream_5kb'].apply(lambda seq:len(seq)==5001)]
numts=numts[numts['downstream_5kb'].apply(lambda seq:len(seq)==5001)]

#create folder for RepeatMasker files
if os.path.isdir('../data/RM_files/')==False:
    os.mkdir('../data/RM_files/')

#function for preparing the output for RepeatMasker analysis, run RepeatMasker itself and process outputs
def repeatmasker(organism_name,output_file):
	try:
		subdf=numts.loc[numts['organism_name']==organism_name]
		#write the sequences into a fasta file
		with open('../data/RM_files/RM_input.fa','w')as outfile:
			subdf.apply(
					lambda row: outfile.write(
							f">{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(len(row['upstream_5kb']))}_{str(int(row['mitochondrial_start']))}\n{row['upstream_5kb']}{row['genomic_sequence'].replace('-','')}{row['downstream_5kb']}\n" #header consists of genomic id+flanking length, genomic and mitochondrial starts (to make sure of uniqueness)
						),
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
#create a df for repeats
#columns are from https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
repeats=pd.DataFrame(
	data=[],
	columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
	)

repeats.to_csv('../results/repeats.csv',index=False)
#make the function work for downstream repeats
organism_names.apply(repeatmasker,args=('../results/repeats.csv',))