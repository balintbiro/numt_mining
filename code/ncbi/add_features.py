#import dependencies
import numpy as np
import pandas as pd

#function for checking mutations
def mutations(seq1,seq2):
    global transversions
    global transitions
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
        transition=[]
        transversion=[]
        no_gap_alignment.loc[:,sites_differ].apply(lambda column: transition.append(1)
                    if (column[0] in purins and column[1] in purins) or (column[0] in pyrimidines and column[1] in pyrimidines)
                    else transversion.append(1),axis=0)
        transversions.append(sum(transversion))
        transitions.append(sum(transition))
    except:
        transversions.append(np.nan)
        transitions.append(np.nan)

#read in numts
numts=pd.read_csv('../data/ncbi_numts_tsne.csv',index_col=0)

#apply the function
transitions=[]
transversions=[]
for index,row in numts.iterrows():
    mutations(row['mitochondrial_sequence'],row['genomic_sequence'])

#add mutations to numts df
numts['transversions']=transversions
numts['transitions']=transitions

#calculate the cumulative sizes of genome parts that contain numts
gpart_cum_sizes=pd.Series(numts['label'].unique()).apply(lambda label: sum(numts.loc[numts['label']==label]['genomic_size'].drop_duplicates()))
gpart_cum_sizes.index=numts['label'].unique()

#calcuate for each numt the ratio of cumulative sizes of genome parts
numt_ratio=numts.apply(lambda row: (row['genomic_length']/gpart_cum_sizes[row['label']])*100,axis=1)

#add numt ratio to numts df
numts['numt_ratio']=numt_ratio

#write out dataframe to csv file
numts.to_csv('../data/ncbi_numts_tsne.csv')