#import dependencies
import os
import numpy as np
import pandas as pd
from subprocess import call
from Bio.SeqUtils import CheckSum
from Bio.SeqUtils import MeltingTemp as mt

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

#get random sequences and numts
random_sequences=pd.read_csv('../data/randoms_rm.csv')
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#function for getting random sequence samples
def get_seq(row):
	sample_start=row['upstream_size']
	if sample_start<0:
		sample_start=0
	sample_end=sample_start+row['sample_size']
	return row['sequence'][sample_start:sample_end]

#get just the sequences
random_sequences['genomic_sequence']=random_sequences.apply(get_seq,axis=1)

#function for calculate GC
def get_GC(sequence):
	if len(sequence)!=0:
		sequence=sequence.upper().replace('-','')
		return (sequence.count('G')+sequence.count('c'))/len(sequence)
	else:
		return np.nan

#calculate GC
random_sequences['GC']=random_sequences['genomic_sequence'].apply(get_GC)

#get flankings of the random sequences
random_sequences['upstream_sequence']=random_sequences.apply(lambda row:row['sequence'][:row['upstream_size']],axis=1)
random_sequences['downstream_start']=random_sequences['upstream_size']+random_sequences['sample_size']
random_sequences['downstream_sequence']=random_sequences.apply(lambda row: row['sequence'][row['downstream_start']:],axis=1)

#add flanking GCs to the random sequences
random_sequences['upstream_GC']=random_sequences['upstream_sequence'].apply(get_GC)
random_sequences['downstream_GC']=random_sequences['downstream_sequence'].apply(get_GC)

#relative positions
random_sequences['rel_start']=random_sequences['genomic_start']/random_sequences['genomic_size']
numts['rel_start']=numts['genomic_start']/numts['genomic_size']

#get infromation entropy
def information_entropy(seq):
	if (type(seq)==str) and (len(seq)!=0):
	    seq=seq.replace('-','').upper()
	    purine,pyrimidine=('A','G'),('C','T')#1,-1 define bases
	    pypy,pypu,pupy,pupu=0,0,0,0#variables for transitions
	    for index in range(0,len(seq)-1):
	        if (seq[index] in pyrimidine) and (seq[index+1] in pyrimidine):
	            pypy+=1
	        elif (seq[index] in pyrimidine) and (seq[index+1] in purine):
	            pypu+=1
	        elif (seq[index] in purine) and (seq[index+1] in purine):
	            pupu+=1
	        else:
	            pupy+=1
	    pyr_trans,pur_trans=(pypy+pypu),(pupu+pypu)
	    try:
	        return -(((pypy/pyr_trans)*np.log(pypy/pyr_trans))
	                +((pypu/pyr_trans)*np.log(pypu/pyr_trans))
	                +((pupu/pur_trans)*np.log(pupu/pur_trans))
	                +((pupy/pur_trans)*np.log(pupy/pur_trans)))
	    except ZeroDivisionError:
	        return np.nan
	else:
		return np.nan

#add information entropies
random_sequences['entropy']=random_sequences['genomic_sequence'].apply(information_entropy)
numts['entropy']=numts['genomic_sequence'].apply(information_entropy)

#flankings information entropies
numts['upstream_entropy']=numts['upstream_5kb'].apply(information_entropy)
numts['downstream_entropy']=numts['downstream_5kb'].apply(information_entropy)
random_sequences['upstream_entropy']=random_sequences['upstream_sequence'].apply(information_entropy)
random_sequences['downstream_entropy']=random_sequences['downstream_sequence'].apply(information_entropy)

#function to add Tms
def melting_temp(seq):
	try:
	    return [
	        mt.Tm_GC(seq,strict=False),#GC Tm
	        mt.Tm_NN(seq,strict=False),#Nearest Neighbor Tm
	        mt.Tm_Wallace(seq,strict=False)
	    ]
	except:
		return 3*[np.nan]

#numt upstream melting temperatures
numt_u_tms=pd.DataFrame(
		data=numts['upstream_5kb'].apply(melting_temp).tolist(),
		columns=['u_TmGC','u_TmNN','u_TmW']
	)
#numts melting temperatures
numt_tms=pd.DataFrame(
		data=numts['genomic_sequence'].apply(melting_temp).tolist(),
		columns=['TmGC','TmNN','TmW']
	)
#numt downstream melting temperatures
numt_d_tms=pd.DataFrame(
		data=numts['downstream_5kb'].apply(melting_temp).tolist(),
		columns=['d_TmGC','d_TmNN','d_TmW']
	)

#upstream random sequence melting temperature
random_u_tms=pd.DataFrame(
		data=random_sequences['upstream_sequence'].apply(melting_temp).tolist(),
		columns=['u_TmGC','u_TmNN','u_TmW']
	)

random_tms=pd.DataFrame(
		data=random_sequences['sequence'].apply(melting_temp).tolist(),
		columns=['TmGC','TmNN','TmW']
	)

random_d_tms=pd.DataFrame(
		data=random_sequences['downstream_sequence'].apply(melting_temp).tolist(),
		columns=['d_TmGC','d_TmNN','d_TmW']
	)

#merge tms to numt df
numts=pd.concat([numts,numt_u_tms,numt_tms,numt_d_tms],axis=1)

#merge tms to randoms df
random_sequences=pd.concat([random_sequences,random_u_tms,random_tms,random_d_tms],axis=1)

#add labels
numts['label']=len(numts)*[1]
random_sequences['label']=len(random_sequences)*[0]

#rename numts GC
numts['GC']=numts['numt_GC']

#read in random repeats
random_repeats=pd.read_csv('../data/random_repeats.csv')

#function for adding repeatranks 
def add_RMranks(row,seq_type):
	try:
		query_name=row['genomic_id']+'_'+str(row.name)
		subdf=random_repeats.loc[random_repeats['query_name']==query_name]
		repeatranks=[]
		if len(subdf)!=0:
			ufil=subdf['piq_begin']<row['upstream_size']
			dfil=subdf['piq_begin']>(row['upstream_size']+row['sample_size'])
			if seq_type=='upstream' and len(subdf[ufil])!=0:
				u_1st_repeat=subdf[ufil]['repeat'].value_counts().index[0]
				u_1st_repeatc=subdf[ufil]['repeat_class'].value_counts().index[0]
				repeatranks+=[u_1st_repeat,u_1st_repeatc]
			elif seq_type=='downstream' and len(subdf[dfil])!=0:
				d_1st_repeat=subdf[dfil]['repeat'].value_counts().index[0]
				d_1st_repeatc=subdf[dfil]['repeat_class'].value_counts().index[0]
				repeatranks+=[d_1st_repeat,d_1st_repeatc]
			else:
				repeatranks+=2*[np.nan]
		else:
			repeatranks+=2*[np.nan]
		return repeatranks
	except:
		return 2*[np.nan]

random_urepeats=pd.DataFrame(data=random_sequences.sample(100).apply(add_RMranks,args=('upstream',),axis=1).tolist(),columns=['u_1st_repeat','u_1st_repeatc'])
random_drepeats=pd.DataFrame(data=random_sequences.sample(100).apply(add_RMranks,args=('downstream',),axis=1).tolist(),columns=['d_1st_repeat','d_1st_repeatc'])

#add random repeats to random sequences
random_sequences=pd.concat([random_sequences,random_urepeats,random_drepeats])

random_sequences[['u_1st_repeat','u_1st_repeatc','d_1st_repeat','d_1st_repeatc']].to_csv('../results/random_rm_ranks.csv')

#select subset dfs
positives=numts[['GC','upstream_GC','downstream_GC','uSW_mean',
       'uSW_median', 'uRMs_count', 'uRMs_lengths', 'dSW_mean', 'dSW_median',
       'dRMs_count', 'dRMs_lengths','rel_start', 'entropy', 'upstream_entropy',
       'downstream_entropy','label','upstream_5kb','genomic_sequence','downstream_5kb',
       'u_TmGC','u_TmNN','u_TmW','TmGC','TmNN','TmW','d_TmGC','d_TmNN','d_TmW']]

negatives=random_sequences[['GC','upstream_GC','downstream_GC','uSW_mean',
       'uSW_median', 'uRMs_count', 'uRMs_lengths', 'dSW_mean', 'dSW_median',
       'dRMs_count', 'dRMs_lengths','rel_start', 'entropy', 'upstream_entropy',
       'downstream_entropy','label','upstream_sequence','genomic_sequence','downstream_sequence',
       'u_TmGC','u_TmNN','u_TmW','TmGC','TmNN','TmW','d_TmGC','d_TmNN','d_TmW']]

#change column names
positives.columns=negatives.columns

features=pd.concat([negatives,positives])

#add sequence sizes
features['upstream_size']=features['upstream_sequence'].apply(lambda seq:len(seq) if type(seq)==str else np.nan)
features['downstream_size']=features['downstream_sequence'].apply(lambda seq:len(seq) if type(seq)==str else np.nan)

features=features[['GC','upstream_GC','downstream_GC','uSW_mean',
       'uSW_median', 'uRMs_count', 'uRMs_lengths', 'dSW_mean', 'dSW_median',
       'dRMs_count', 'dRMs_lengths','rel_start', 'entropy', 'upstream_entropy',
       'downstream_entropy','label','upstream_size','downstream_size',
       'u_TmGC','u_TmNN','u_TmW','TmGC','TmNN','TmW','d_TmGC','d_TmNN','d_TmW']]
features.to_csv('../data/ml_features.csv',index=False)