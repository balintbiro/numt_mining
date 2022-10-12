#import dependencies
import os
import numpy as np
import pandas as pd
from subprocess import call
import matplotlib.pyplot as plt

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

repeats=pd.read_csv('../results/longest_chr_RM.csv')

#filtering numts df - to get the required species and chrs
numts=numts.loc[numts['organism_name'].isin(organism_names)]#required species
numts=numts.loc[numts['genomic_id'].isin(longest_chr_ids)]#required chr ids

#get query ids - in that way it is possible to synchronize numts and repeats
numts['query_name']=numts.apply(
	lambda row: f"{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(int(row['mitochondrial_start']))}", #header consists of genomic id+flanking length, genomic and mitochondrial starts (to make sure of uniqueness)
	axis=1
)

#get relelvant information for sliding window
numts=numts[['query_name','upstream_5kb','genomic_length','downstream_5kb']]
numts['upstream_length']=numts.apply(lambda row: len(row['upstream_5kb']),axis=1)
numts['downstream_length']=numts.apply(lambda row: len(row['downstream_5kb']),axis=1)
numts=numts[['query_name','upstream_length','genomic_length','downstream_length']]

#get relevant repeat rows
repeats=repeats[['query_name','piq_begin','piq_end','repeat','repeat_class']]

#get numt ranges
numt_ranges=numts.apply(lambda row: np.arange(5001,5001+row['genomic_length']),axis=1)
numt_ranges.index=numts['query_name']

#function for getting ranges of repeat classes
def get_rm_range(row,repeat_class):
	global repeat_ranges
	query_id=row['query_name']
	repeat_subdf=repeats.loc[repeats['query_name']==query_id]
	if repeat_class in repeat_subdf['repeat_class'].tolist():
		repeat_subdf=repeat_subdf.loc[repeat_subdf['repeat_class']==repeat_class]
		repeat_range=np.concatenate(repeat_subdf.apply(lambda row:np.arange(row['piq_begin'],row['piq_end']),axis=1).tolist()).tolist()
		repeat_ranges.append([query_id,repeat_class,repeat_range])

#global variable
repeat_ranges=[]

#get ranges for every repeat class
for repeat_class in repeats['repeat_class'].unique():
	numts.apply(get_rm_range,args=(repeat_class,),axis=1)

#create df from repeat ranges
repeat_ranges=pd.DataFrame(repeat_ranges)
repeat_ranges.columns=['query_name','repeat_class','repeat_range']
repeat_ranges=repeat_ranges.set_index('query_name',drop=True)

#function for sliding window
def sliding_window(query_name,window_size):
	global frequencies
	repeat_classes=repeat_ranges.loc[query_name]['repeat_class'].tolist()
	for repeat_class in repeat_classes:
		rm_frequency=[]
		numt_frequency=[]
		numt_range=set(numt_ranges[query_name])
		repeat_range=repeat_ranges.loc[query_name]
		repeat_range=set(repeat_range.loc[repeat_range['repeat_class']==repeat_class]['repeat_range'].tolist()[0])
		for step in np.linspace(0,5000+list(numt_range)[-1],2500,dtype=int):
			window=set(np.arange(step,step+window_size))
			rm_frequency.append(len(repeat_range&window))
			numt_frequency.append(len(numt_range&window))
		frequencies+=[[query_name,repeat_class,'rm',rm_frequency]]
		if [query_name,repeat_class,'numt',numt_frequency] not in frequencies:
			frequencies+=[[query_name,'numt','numt',numt_frequency]]

#global variable
frequencies=[]

pd.Series(numt_ranges.index).apply(sliding_window,args=(50,))

#load sliding windows into a df
sliding_windows=pd.DataFrame(frequencies)
sliding_windows.columns=['query_name','repeat_class','type','range']
sliding_windows=sliding_windows.drop_duplicates(subset=['query_name','repeat_class'],keep='first')
sliding_windows['organism_name']=sliding_windows['query_name'].apply(lambda query_name:longest_chr_ids[longest_chr_ids=='_'.join(query_name.split('_')[:2])].index[0])
sliding_windows.to_csv('../results/sliding_windows.csv')

#create directory for sliding window figures
if os.path.isdir('../results/sliding_windows/')==False:
	os.mkdir('../results/sliding_windows/')

#function for plotting the results
def plotter(query_name):
    subdf=sliding_windows.loc[sliding_windows['query_name']==query_name]
    organism_name=subdf['organism_name'].tolist()[0]
    fig,axs=plt.subplots(1,1,figsize=(10,6))
    for index, row in subdf.iterrows():
    	if row['repeat_class']=='numt':
    		axs.plot(row['range'],label=row['repeat_class'],lw=3,alpha=.3)
    	else:
        	axs.plot(row['range'],label=row['repeat_class'])
    plt.legend()
    plt.savefig(f'../results/sliding_windows/{organism_name}_{query_name}.png',dpi=100)

pd.Series(sliding_windows['query_name'].unique()).apply(plotter)
