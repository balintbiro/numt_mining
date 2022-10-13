#importing dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#reading basic data
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#select flankings
flankings=numts[['upstream_5kb','genomic_sequence','downstream_5kb']]

#filter rows where flankings have the appropriate lengths
upstream_fil=flankings['upstream_5kb'].apply(lambda upstream:type(upstream)==str and len(upstream)>200)
downstream_fil=flankings['downstream_5kb'].apply(lambda downstream:type(downstream)==str and len(downstream)>200)
flankings=flankings[upstream_fil][downstream_fil]

#define function for sequence creation
def get_seq(row):
	global sequences
	upstream_part=row['upstream_5kb'][-200:].upper()
	numt_part=row['genomic_sequence'].replace('-','').upper()[:10]+row['genomic_sequence'].replace('-','').upper()[-10:]
	downstream_part=row['downstream_5kb'][:200].upper()
	sequences.append(list(upstream_part+numt_part+downstream_part))

#global variable
sequences=[]

#make the function work
flankings.apply(get_seq,axis=1)

#create dataframe from sequences
positions=pd.DataFrame(sequences)

#function for counting nucleotides in each position
def get_count(column,nucleotide):
    return ''.join(column.tolist()).count(nucleotide)

#get counts for each positions for each nucleotides
A_counts=positions.apply(get_count,args=('A',),axis=0)
C_counts=positions.apply(get_count,args=('C',),axis=0)
G_counts=positions.apply(get_count,args=('G',),axis=0)
T_counts=positions.apply(get_count,args=('T',),axis=0)

#visualize the result
plt.style.use('fivethirtyeight')
fig,axs=plt.subplots(1,1,figsize=(10,8))
x=np.linspace(-200,200,len(A_counts),dtype=int)
axs.plot(x,A_counts,'blue',linewidth=3,label='A')
axs.plot(x,C_counts,'black',linewidth=3,label='C')
axs.plot(x,G_counts,'pink',linewidth=3,label='G')
axs.plot(x,T_counts,'red',linewidth=3,label='T')
plt.fill_between([-10,10],axs.get_ylim()[0],axs.get_ylim()[1],alpha=.5,color='lightblue')
plt.xticks([-200,-150,-100,-50,0,50,100,150,200])
axs.set_xticklabels([-200,-150,-100,-50,'numt',50,100,150,200])
plt.legend()
plt.tight_layout()
plt.savefig('../results/base_composition.png',dpi=400)