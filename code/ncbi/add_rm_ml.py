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
organism_names.apply(repeatmasker,args=('../data/random_repeats.csv',))

#get sequences and repeats
repeats=pd.read_csv('../data/random_repeats.csv')
sequences=pd.read_csv('../data/ml_input_randoms.csv')

#get indices
indices=repeats['query_name'].apply(lambda query_name: int(query_name.split('_')[-1]))

#add upstream, downstream sizes and sample sizes (which referes to numt sizes)-use the query names from repeats!
repeats['upstream_size']=sequences.iloc[indices]['upstream_size'].tolist()
repeats['downstream_size']=sequences.iloc[indices]['upstream_size'].tolist()
repeats['sample_size']=sequences.iloc[indices]['sample_size'].tolist()

#select upstream and downstream repeats
upstream_repeats=repeats[repeats['piq_begin']<repeats['upstream_size']]
downstream_repeats=repeats[repeats['piq_begin']>(repeats['upstream_size']+repeats['sample_size'])]

#add repeats to sequences df
#function for adding different features of RepeatMasker to the sequence df
def add_RM(query_name,sequence_type):
	if sequence_type=='upstream':
		flanking_repeats=upstream_repeats
	else:
		flanking_repeats=downstream_repeats
	try:
		subdf=flanking_repeats.loc[flanking_repeats['query_name']==query_name]
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

#get query names
query_names=sequences.apply(lambda row: row['genomic_id']+'_'+str(row.name),axis=1)

#apply the function to query names
uRM_features=pd.DataFrame(query_names.apply(add_RM,args=('upstream',)).tolist())
uRM_features.columns=['uSW_mean','uSW_median','uRMs_count','uRMs_lengths']

dRM_features=pd.DataFrame(query_names.apply(add_RM,args=('downstream',)).tolist())
dRM_features.columns=['dSW_mean','dSW_median','dRMs_count','dRMs_lengths']

#add repeats to sequences
output_df=pd.concat([sequences,uRM_features,dRM_features],axis=1)

#write updated repeats 
output_df.to_csv('../data/randoms_rm.csv',index=False)