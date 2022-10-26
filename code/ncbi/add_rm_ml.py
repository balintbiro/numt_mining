#import dependencies
import os
import pandas as pd
from subprocess import call

#reading in randomly selected sequences
randoms=pd.read_csv('../data/ml_input_randoms.csv')

#create folder for RepeatMasker files
if os.path.isdir('../data/RM_files/')==False:
    os.mkdir('../data/RM_files/')

#function for preparing the output for RepeatMasker analysis, run RepeatMasker itself and process outputs
def repeatmasker(organism_name,output_file):
	try:
		subdf=randoms.loc[randoms['organism_name']==organism_name]
		#write the sequences into a fasta file
		with open('../data/RM_files/RM_input.fa','w')as outfile:
			subdf.apply(
					lambda row: outfile.write(
							f">{row['genomic_id']}_{str(row.name)}\n{row['sequence']}\n" #header consists of genomic id+frow index (to make sure of uniqueness)
						)
					if type(row['sequence'])==str #in some cases there is no flanking region and so np.nans are present
					else print(f"{organism_name} doesnt have sequence under the genomic ID of {row['genomic_id']} in the position of {row.name}"),
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
randoms['organism_name']=randoms['organism_name'].apply(lambda organism_name: organism_name.lower())
organism_names=randoms['organism_name'].drop_duplicates()

#create an empty df for random repeats
random_repeats=pd.DataFrame(
		data=[],
		columns=[
			'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
		]
	)

#write empty df
random_repeats.to_csv('../data/random_repeats.csv',index=False)

#make the function work
organism_names[:3].apply(repeatmasker,args=('../data/random_repeats.csv',))