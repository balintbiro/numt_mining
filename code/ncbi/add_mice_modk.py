#import dependencies
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

#get filenames
filenames=pd.Series(os.listdir('../data/mice_strain_numts/')).apply(
    lambda name: name if name[0]!='.'
    else np.nan
)

#dropna
filenames=filenames.dropna()

#get indices
indices=filenames.apply(lambda name:'_'.join(name.split('_')[:3]) if '.' not in name.split('_')[2]
               else '_'.join(name.split('_')[:2]))

#read in dataframes
dfs=filenames.apply(
    lambda filename: pd.read_csv(f'../data/mice_strain_numts/{filename}')
)
dfs.index=indices

#add labels
def add_labels(df):
    global indexer
    df['name']=len(df)*[indices.tolist()[indexer]]
    df['label']=len(df)*[indexer]
    indexer+=1
    
indexer=0
dfs.apply(add_labels)

#merge dfs into one merged df
numts=pd.concat(dfs.tolist())

#function for calculating modified Kimura2 parameter
#https://link.springer.com/article/10.1007/s00239-018-9885-1
def modK2(seq1,seq2):
    global transversions
    global transitions
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
    transition=[]
    transversion=[]
    Ks=0
    try:
        no_gap_alignment.loc[:,sites_differ].apply(lambda column: transition.append(1)
                if (column[0] in purins and column[1] in purins) or (column[0] in pyrimidines and column[1] in pyrimidines)
                else transversion.append(1),axis=0)
        S=len(no_gap_alignment.loc[:,~sites_differ].columns.values)/len(alignment.columns.values)
        P=sum(transition)/len(alignment.columns.values)
        Q=sum(transversion)/len(alignment.columns.values)
        K=(3/4)*w*np.log(w)-(w/2)*np.log(S-P)*np.sqrt(S+P-Q)
        Ks=K
        transversions.append(sum(transversion))
        transitions.append(sum(transition))
    except:
        Ks=Ks
        transversions.append(0)
        transitions.append(0)
    return Ks

#calculate Kimura2 parameter
Kimura2s=[]
transitions=[]
transversions=[]
for index,row in numts.iterrows():
    Kimura2s.append(modK2(row['mt_sequence'],row['g_sequence']))

#add Kimura2 Distances to numts dataframe
numts['modk2']=Kimura2s
numts['transitions']=transitions
numts['transversions']=transversions

#calculate GC content of numts and add them to numts df
numts['GC']=numts['g_sequence'].apply(
    lambda seq: (seq.upper().count('G')+seq.upper().count('C'))/len(seq.replace('-',''))
    )

#export df to csv
numts.to_csv('../data/mice_numts.csv')