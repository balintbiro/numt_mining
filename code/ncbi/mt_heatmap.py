#import dependencies
import numpy as np
import pandas as pd
from Bio import AlignIO

#read in numts
numts=pd.read_csv('../data/ncbi_numts_p26.csv')
numts=numts[['organism_name','mitochondrial_start','mitochondrial_length']]

#calculate numt ranges
numts['numt_ranges']=numts.apply(
		lambda row: list(np.arange(row['mitochondrial_start'],row['mitochondrial_start']+row['mitochondrial_length'])),
		axis=1
	)

#get organism names
organism_names=pd.Series(numts['organism_name'].unique())

#function for concatenating numt ranges for each organism
def get_numt_ranges(organism_name):
	subdf=numts.loc[numts['organism_name']==organism_name]
	return sum(subdf['numt_ranges'].tolist(),[])

#make the function work
numt_ranges=organism_names.apply(get_numt_ranges)
numt_ranges.index=organism_names

#read in alignments
msa=AlignIO.read(
		open('../results/aligned_mtDNAs.fa'),
		format='fasta'
	)

#add ids to msa
msa_ids=[]
for record in msa:
	msa_ids.append(record.id.lower())
msa=pd.DataFrame(msa)
msa.index=msa_ids

#transform msa-do not show individual nucleotides but their position in mitochondrion
def transform_msa(organism_name):
	alignment_range=msa.loc[organism_name]
	alignment_range[alignment_range!='-']=np.arange(0,len(alignment_range[alignment_range!='-']))
	msa.loc[organism_name]=alignment_range

#function for getting numt_coverages
def get_numt_coverage(organism_name,output_file):
	numt_range=pd.Series(numt_ranges[organism_name]).value_counts()
	numt_range['-']=-1#add np.nan for gaps
	msa_range=msa.loc[organism_name]
	numt_coverage=list(numt_range.reindex(msa_range.values).fillna(0).values)
	numt_coverage=[organism_name]+numt_coverage
	numt_coverage=pd.DataFrame([numt_coverage])
	numt_coverage.to_csv(output_file,mode='a',index=False,header=False)

#transform msa from alignment to mt positions
organism_names.apply(transform_msa)

#create empty df
heatmap=pd.DataFrame(data=[],columns=np.arange(0,msa.shape[1]))
#write it to csv
heatmap.to_csv('../results/heatmap.csv',index=False)

#transform df to numt coverages
organism_names.apply(get_numt_coverage,args=('../results/heatmap.csv',))

####################################################################################################
#                                       visualisation                                              #
#                        for just the visualization, copy the code from here                       #
####################################################################################################

#read in file
heatmap_input=pd.read_csv('../results/heatmap.csv')

#visualize results
mask=heatmap_input<0

#create color codes
taxonomy_data=pd.read_csv('../data/taxonomy_data.csv',index_col=0)
taxonomy_data=taxonomy_data.drop_duplicates(subset='organism_name')
taxonomy_data=taxonomy_data.set_index('organism_name',drop=True)
order_dict=taxonomy_data['order'].dropna()
orders=order_dict.reindex(heatmap_input.index)
colors=pd.Series(     
    data=['#1f77b4ff','#aec7e8ff','#ff7f0eff','#ffbb78ff','#2ca02cff','#98df8aff','#d62728ff','#ff9896ff',          
          '#9467bdff','#c5b0d5ff','#8c564bff','#c49c94ff','#e377c2ff','#f7b6d2ff','#7f7f7fff','#c7c7c7ff',          
          '#bcbd22ff','#dbdb8dff','#FFFFFF','#FFFFFF','#FFFFFF'],     
    index=['Carnivora','Primates','Chiroptera','Artiodactyla','Rodentia','Pilosa','Eulipotyphla','Cingulata',            
           'Microbiotheria','Macroscelidea','Dermoptera','Didelphimorphia','Proboscidea','Pholidota','Lagomorpha',           
           'Monotremata','Diprotodontia','Dasyuromorphia','Perissodactyla','Tubulidentata',np.nan] 
)
order_colors=colors.reindex(orders.values).tolist()

#robust means that no outliers are displayed
heatmap=sns.clustermap(heatmap_input,cmap='rocket_r',mask=mask,robust=True,row_colors=order_colors,cbar_pos=(-.1, .2, .03, .4))
ax=heatmap.ax_heatmap
cbar=ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
ax.set_xticks([])
ax.set_yticks([])
heatmap.savefig('../results/heatmap.png',dpi=400)
