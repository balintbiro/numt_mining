#import depepndencies
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#load dataframe
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

########################################################################
#numt sizes
#calculate relative (to mt size) numt size
rel_numt_size=numts['genomic_length']/(numts['mtDNA_size (Mb)']*1000000)

#visualize reuslts
fig,ax=plt.subplots(1,1)
ax.hist(rel_numt_size[rel_numt_size<0.1],color='brown',alpha=0.35,bins=100)
ax.set_xlabel('Relative numt size')
ax.set_ylabel('numt count')

l,b,w,h = .65, .35, .2, .45#left, bottom, width and height
ax1=fig.add_axes([l,b,w,h])
boxprops = dict(linestyle='-', linewidth=2, color='black')
ax1.boxplot(numts['genomic_length'],showfliers=False,boxprops=boxprops)
ax1.spines['top'].set_color('0.5')
ax1.spines['bottom'].set_color('0.5')
ax1.spines['left'].set_color('0.5')
ax1.spines['right'].set_color('0.5')
ax1.set_xticklabels([])
ax1.set_ylabel('numt size (bp)')
ax1.set_title('Abs. numt size')
plt.tight_layout()
plt.savefig('../results/numt_lengths.png',dpi=400)
########################################################################

#numt content correlations
numt_per_gid=numts['genomic_id'].drop_duplicates().apply(
		lambda genomic_id: len(numts.loc[numts['genomic_id']==genomic_id])
	)
gid_sizes=numts.drop_duplicates(subset='genomic_id')['genomic_size']/1000000
numtbp_per_gid=numts['genomic_id'].drop_duplicates().apply(
		lambda genomic_id: sum(numts.loc[numts['genomic_id']==genomic_id]['genomic_length'])
	)/1000

#visualize
fig,axs=plt.subplots(2,1,figsize=(7,6),sharex=True)
axs[0].plot(gid_sizes,numt_per_gid,'o',alpha=.15)
axs[0].set_ylabel('Number of numts')
a0,b0=np.polyfit(gid_sizes,numt_per_gid,1)
axs[0].plot(gid_sizes,a0*gid_sizes+b0,color='red',lw=1)

axs[1].plot(gid_sizes,numtbp_per_gid,'o',color='green',alpha=.15)
a1,b1=np.polyfit(gid_sizes,numtbp_per_gid,1)
axs[1].plot(gid_sizes,a1*gid_sizes+b1,color='red',lw=1)
axs[1].set_xlabel('Size (Mb)')
axs[1].set_ylabel('numt bps (kb)')
plt.tight_layout()
plt.savefig('../results/sizes_corr.png',dpi=400)
########################################################################

#relative GCs
GC_input=numts.dropna(subset=['genomic_sequence','mitochondrial_sequence','gDNA_GC','mtDNA_GC'])
#calculate numt GC and normalize it with the whole genome GC
gseq_GC=GC_input['genomic_sequence'].apply(lambda seq: (seq.upper().replace('-','').count('G')+seq.upper().replace('-','').count('C'))/len(seq.replace('-','')))
grel_GC=(gseq_GC)*100/GC_input['gDNA_GC']

#calculate numt corresponding mt seq GC and normalize it with whole mt GC
mtseq_GC=GC_input['mitochondrial_sequence'].apply(lambda seq: (seq.upper().replace('-','').count('G')+seq.upper().replace('-','').count('C'))/len(seq.replace('-','')))
mtrel_GC=(mtseq_GC)*100/GC_input['mtDNA_GC']

#visualize results
fig,axs=plt.subplots(1,1)
axs.hist(grel_GC,bins=500,color='blue',alpha=.3,label='grelGC')
axs.hist(mtrel_GC,bins=500,color='green',alpha=.3,label='mtrelGC')
plt.legend()
plt.xlabel('Relative GC content')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig('../results/relGCs.png',dpi=400)
