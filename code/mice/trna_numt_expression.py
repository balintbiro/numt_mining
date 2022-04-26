#import required modules
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from subprocess import call

#reading mm numts
numts=pd.read_csv('..results/numts_with_same_mt/Mus_musculus_numts.csv')

#reading mmm mt annotation
mt_annot=pd.read_csv('/../data/mt_annotations/mm.bed',
                    sep='\t',header=None)
mt_annot.columns=['MT','start','end','name','sig','strand']
mt_annot['length']=mt_annot.apply(lambda row: row['end']-row['start'],axis=1)

#get trna annotations
mt_trnas=mt_annot.loc[mt_annot['name'].str.contains('trn')]
trna_ranges=mt_trnas.apply(lambda row: np.arange(row['start'],row['end']),axis=1)
trna_ranges.index=mt_annot.loc[mt_annot['name'].str.contains('trn')]['name']

#get numt ranges
numt_ranges=numts.sort_values(by='mt_start').apply(lambda row: np.arange(row['mt_start'],row['mt_start']+row['mt_length']),axis=1)

#define function to check the intersections of numts and trnas
def tRNA_numt(tRNA_range):
    trial=[]
    for item in numt_ranges:
        trial.append(len(set(tRNA_range)&set(item)))
    return pd.Series(trial)
    
intersections=trna_ranges.apply(tRNA_numt)

#replace 0-s to nans
intersections.columns=numt_ranges.index.values
intersections.replace(to_replace=0, value=np.nan, inplace=True)

#add trna length
intersections=intersections.dropna(how='all',axis=1)
intersections['trna_length']=trna_ranges.apply(lambda trna_range:len(trna_range))

#drop artefact
intersections=intersections.drop(labels=[0],axis=1)

#get the tRNA pseudogenes (aka numts)
trna_numts=intersections.apply(lambda row: row.tolist().index(row['trna_length']),axis=1).apply(lambda index: numts.loc[index]['g_sequence'].upper())

#get mm mt seq
mt_seq=''
for seq_record in SeqIO.parse('../data/Mus_musculus/mt.fa', "fasta"):
    mt_seq=str(seq_record.seq)

#get tRNA sequences
trna_seqs=mt_trnas.apply(lambda row: mt_seq[row['start']:row['end']],axis=1)
trna_seqs.index=mt_trnas['name']

if os.path.isdir('../data/mice_trna_sequences/')==False:
	os.mkdir('../data/mice_trna_sequences/')

if os.path.isdir('../data/mice_trna_numts/')==False:
	os.mkdir('../data/mice_trna_numts/')

#function for writing output files
def write_output(sequence,sequences,output_dir):
	global indexer
	filename=f'{sequences.index.values[indexer]}.fa'
	with open(os.path.join(output_dir+filename),'w')as outfile:
		outfile.write('>'+sequences.index.values[indexer]+'\n'+sequence)

indexer=0
trna_numts.apply(write_output,args((trna_numts,'../data/mice_trna_numts/',)))

indexer=0
