#importing required modules
import os
import numpy as np
import pandas as pd
from Bio import AlignIO
from subprocess import call

#function for performing alignment with EMBOSS needle
def alignment(filename):
	if os.path.isdir('../data/mt_pairwise_alignments/')==False:
		os.mkdir('../data/mt_pairwise_alignments')
	#outfilename=f'{filename.split('.')[0].split('_')[-1]}.afa'
	call(f'needle -asequence ../data/mice_mts/Mus_musculus.fa -bsequence ../data/mice_mts/{filename} -gapopen 10 -gapextend 0.5 -outfile ../data/mt_pairwise_alignments/{filename}',shell=True)

#execute alignment
files=pd.Series(os.listdir('../data/mice_mts/'))
files.index=files
files=files.drop('Mus_musculus.fa')
#files.apply(alignment)

if os.path.isdir('../results/modK2_parameters/')==False:
	os.mkdir('../results/modK2_parameters/')

#function for calculating modified Kimura2 paramter based on:
#https://link.springer.com/article/10.1007/s00239-018-9885-1
def modK2_dist(filename,window_size,step_size):
    alignment=AlignIO.read(f'../data/mt_pairwise_alignments/{filename}','emboss')
    alignment=pd.DataFrame(alignment)
    Ks=[]
    for i in np.arange(0,len(alignment.columns.values),step_size):
        if i+window_size<len(alignment.columns.values):
            sub_alignment=alignment.loc[:,np.arange(i,i+window_size)]
            #get number of nucleotides and probability of nucleotides
            number_of_nucleotides=(len(''.join(sub_alignment.loc[0].tolist()).replace('-',''))+
                                   len(''.join(sub_alignment.loc[1].tolist()).replace('-','')))
            w=number_of_nucleotides/(len(sub_alignment.columns.values)*2)#probability of nucleotides (not gaps!)
            #delete gaps from alignment
            sub_alignment=sub_alignment.apply(lambda row: row.replace('-',np.nan),axis=1)
            no_gap_alignment=sub_alignment.dropna(axis=1)
            #get sites where nucleotides are not the same
            sites_differ=no_gap_alignment.apply(lambda column: column[0]!=column[1],axis=0)
            purins=['A','G']
            pyrimidines=['C','T']
            transitions=[]
            transversions=[]
            try:
                no_gap_alignment.loc[:,sites_differ].apply(lambda column: transitions.append(1)
                        if (column[0] in purins and column[1] in purins) or (column[0] in pyrimidines and column[1] in pyrimidines)
                        else transversions.append(1),axis=0)
                S=len(no_gap_alignment.loc[:,~sites_differ].columns.values)/len(sub_alignment.columns.values)
                P=sum(transitions)/len(sub_alignment.columns.values)
                Q=sum(transversions)/len(sub_alignment.columns.values)
                K=(3/4)*w*np.log(w)-(w/2)*np.log(S-P)*np.sqrt(S+P-Q)
                Ks.append(K)
            except:
                Ks.append(0)
    plt.style.use('fivethirtyeight')
    fig,axs=plt.subplots(1,1)
    axs.plot(Ks)
    axs.set_xlabel('Mitochondrial_nucleotides (bp)')
    axs.set_ylabel('mod K2 distance')
    plt.tight_layout()
    plt.savefig(f"../results/modK2_parameters/{filename.split('.')[0].split('_')[-1]}_modK2.png",dpi=200)

#get modK2 parameters
pd.Series(os.listdir('../data/mt_pairwise_alignments/')).apply(modK2_dist, args=(1000,10,))