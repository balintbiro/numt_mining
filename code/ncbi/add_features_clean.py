#import dependencies
import os
import  numpy as np
import pandas as pd
from subprocess import call

#reading different dfs
dfs=pd.Series(os.listdir('../data/alignments/')).apply(lambda csv_file: pd.read_csv(f'../data/alignments/{csv_file}'))

#get organism names
organism_names=pd.Series(os.listdir('../data/alignments/'))
organism_names=organism_names.apply(lambda name: '_'.join(name.split('_')[:2]).lower()).tolist()

#add organism_names to dfs
for index, organism_name in enumerate(organism_names):
	dfs[index]['organism_name']=len(dfs[index])*[organism_name]

#merging different dfs into a common df
numts=pd.concat(dfs.tolist())

#write common df into a csv file
numts.to_csv('../data/ncbi_numts.csv')

#function for calculating modified Kimura2 parameter
#https://link.springer.com/article/10.1007/s00239-018-9885-1
def modK2(seq1,seq2):
    try:
        seq1=seq1.upper()
        seq2=seq2.upper()
        alignment=pd.DataFrame(
            [
                list(seq1),
                list(seq2)
            ]
        )
        #get number of nucleotides and probability of nucleotides
        number_of_nucleotides=(len(''.join(alignment.loc[0].tolist()).replace('-',''))+
                               len(''.join(alignment.loc[1].tolist()).replace('-','')))
        #probability of nucleotides not gaps
        w=number_of_nucleotides/(len(alignment.columns.values)*2)
        #delete gaps from alignment
        alignment=alignment.apply(lambda row: row.replace('-',np.nan),axis=1)
        no_gap_alignment=alignment.dropna(axis=1)
        #get sites where nucleotides are not the same
        sites_differ=no_gap_alignment.apply(lambda column: column[0]!=column[1],axis=0)
        purins=['A','G']
        pyrimidines=['C','T']
        transitions=[]
        transversions=[]
        Ks=0
        try:
            no_gap_alignment.loc[:,sites_differ].apply(lambda column: transitions.append(1)
                    if (column[0] in purins and column[1] in purins) or (column[0] in pyrimidines and column[1] in pyrimidines)
                    else transversions.append(1),axis=0)
            S=len(no_gap_alignment.loc[:,~sites_differ].columns.values)/len(alignment.columns.values)
            P=sum(transitions)/len(alignment.columns.values)
            Q=sum(transversions)/len(alignment.columns.values)
            K=(3/4)*w*np.log(w)-(w/2)*np.log(S-P)*np.sqrt(S+P-Q)
            Ks=K
        except:
            Ks=Ks
        return Ks
    except:
        return np.nan

#calculate Kimura2 parameter
Kimura2s=[]
for index,row in numts.iterrows():
    Kimura2s.append(modK2(row['mitochondrial_sequence'],row['genomic_sequence']))

#add Kimura2 Distances to numts dataframe
numts['modk2']=Kimura2s

#function for calculating GC content
def get_GC(seq):
	try:
		return (seq.upper().count('G')+seq.upper().count('C'))/len(seq.replace('-',''))
	except:
		return np.nan

#calculate GCs and add them to the merged df
numts['numt_GC']=numts['genomic_sequence'].apply(get_GC)
numts['upstream_GC']=numts['upstream_5kb'].apply(get_GC)
numts['downstream_GC']=numts['downstream_5kb'].apply(get_GC)

#export csv p4 means that 4 patterns have been added (Kimura2, numtGC, flankings GC)
numts.to_csv('../data/ncbi_numts_p4.csv',index=False)