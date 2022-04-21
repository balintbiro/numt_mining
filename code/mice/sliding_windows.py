#import required modules
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

#load repeatmasker df
rmsk=pd.read_csv('../data/mm_rmsk', sep='\t')
rmsk=rmsk.set_index('chrom')
rmsk['repclass']=rmsk['name'].apply(lambda rm_name:rm_name.split('#')[1])

#load numts
numts=pd.read_csv('../data/Mus_musculus_numts.csv')
numts=numts.set_index('g_id')

#get just the crhomosomes
rmsk=rmsk.loc[['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19']]
numts=numts.loc[['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '2',
       '3', '4', '5', '6', '7', '8', '9']]
rmsk=rmsk.loc[['chr18']]

#cut chr off from chr ids
rmsk.index=pd.Series(rmsk.index.values).apply(lambda chr_id:chr_id.split('chr')[1])

numt_test=numts.loc['18'].sort_values(by='g_start')
numt_array=set(np.concatenate(numt_test.apply(lambda row: np.arange(row['g_start'],(row['g_start']+row['g_length'])),axis=1)))

#function for sliding window
def sliding_window(repname,chr_size,window_size,step_size,numt_array):
    plt.style.use('fivethirtyeight')
    rmsk_test=rmsk.loc[rmsk['name']==repname]
    rmsk_array=set(np.concatenate(rmsk_test.apply(lambda row: np.arange(row['chromStart'],row['chromEnd']),axis=1)))
    repeats=[]
    numts=[]
    for i in np.arange(0,chr_size,step_size):
        window=set(np.arange(i,i+window_size))
        repeats.append((len(window&rmsk_array)/window_size)*100)
        numts.append((len(window&numt_array)/window_size)*100)
    if pearsonr(repeats,numts)[0]<-0.1 or pearsonr(repeats,numts)[0]>0.1:
        fig,axs=plt.subplots(2,1)
        axs[0].plot(numts,linewidth=0.8)
        axs[1].plot(repeats,linewidth=0.8)
        axs[0].fill(np.arange(0,len(numts)),numts,alpha=0.3)
        axs[1].fill(np.arange(0,len(repeats)),repeats,alpha=0.3)
        axs[0].set_title(repname+' '+str(round(pearsonr(repeats,numts)[0],4)))
        repname=re.sub(r'[\\/*?:"<>|]',"",repname)
        plt.savefig(f'../data/sliding_windows/{repname}.png')
    else:
        print(repname)

pd.Series(rmsk['name'].unique()).apply(sliding_window,args=(numt_test['g_size'][0],50000,50000,numt_array))