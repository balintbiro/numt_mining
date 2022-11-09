#import dependencies
import os
import pandas as pd
from subprocess import call
import matplotlib.pyplot as plt

#reading in merged numts df
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#filter df
numts=numts.dropna(subset=['upstream_5kb','downstream_5kb'])
numts=numts[['organism_name','genomic_id','genomic_start','upstream_5kb','genomic_length','downstream_5kb','mitochondrial_start','genomic_sequence']]
numts=numts[numts['upstream_5kb'].apply(lambda seq:len(seq)==5001)]
numts=numts[numts['downstream_5kb'].apply(lambda seq:len(seq)==5001)]

#create folder for RepeatMasker files
if os.path.isdir('../data/RM_files/')==False:
    os.mkdir('../data/RM_files/')

#function for preparing the output for RepeatMasker analysis, run RepeatMasker itself and process outputs
def repeatmasker(organism_name,output_file):
	try:
		subdf=numts.loc[numts['organism_name']==organism_name]
		#write the sequences into a fasta file
		with open('../data/RM_files/RM_input.fa','w')as outfile:
			subdf.apply(
					lambda row: outfile.write(
							f">{row['genomic_id']}_{str(int(row['genomic_start']))}_{str(int(row['genomic_length']))}_{str(len(row['upstream_5kb']))}_{str(int(row['mitochondrial_start']))}\n{row['upstream_5kb']}{row['genomic_sequence'].replace('-','')}{row['downstream_5kb']}\n" #header consists of genomic id+flanking length, genomic and mitochondrial starts (to make sure of uniqueness)
						),
					axis=1
				)
		#change directory and call RepeatMasker
		call(f"""cd ../data/RM_files; RepeatMasker -species {organism_name.replace('_',' ')} RM_input.fa""",shell=True)
		#transform RepeatMasker.out file into a tab separated file for better handling
		call("""cat ../data/RM_files/RM_input.fa.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' > ../data/RM_files/RM_input.tab""",shell=True)
		#add columns https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
		repeats=pd.read_csv('../data/RM_files/RM_input.tab',sep='\t',skiprows=2,header=None)
		#add organism name to repeats df
		repeats['organism_name']=len(repeats)*[organism_name]
		#add current repeats to an existing output file
		repeats.to_csv(output_file,mode='a',index=False,header=False)
		#clear the directory
		call('rm -r ../data/RM_files/*',shell=True)
	except:
		pass

#get organism names
organism_names=pd.Series(os.listdir('../data/alignments/')).apply(lambda name: '_'.join(name.split('_')[:2]).lower())

######################################################################################################
#create a df for repeats
#columns are from https://www.repeatmasker.org/webrepeatmaskerhelp.html How to read the results
repeats=pd.DataFrame(
	data=[],
	columns=[
	'SW_score','div_perc','del_perc','ins_perc','query_name','piq_begin','piq_end','piq_left','matching','repeat',
	'repeat_class','pir_begin','pir_end','pin_left','ID','higher_score_match','organism_name'
	]#piq means Position In Query while pir means Position In Repeat
	)

repeats.to_csv('../results/repeats.csv',index=False)
#make the function work for downstream repeats
organism_names.apply(repeatmasker,args=('../results/repeats.csv',))

#read in repeats
repeats=pd.read_csv('../results/repeats.csv')

#create column for numt lengths
repeats['numt_lengths']=repeats['query_name'].apply(lambda query_name: int(query_name.split('_')[-3]))

#get the 40 downstream nucleotides of a given numt plus the 5kb flanking
repeats['downstream_pos']=repeats.apply(lambda row: (5000+(row['numt_lengths']-40)),axis=1)

#filter downstream repeats - it is different for every numt because of the size variability
downstream_repeats=repeats[repeats['piq_begin']>repeats['downstream_pos']]
downstream_repeats['dpiq_begin']=downstream_repeats['piq_begin']-downstream_repeats['numt_lengths']

#function for plotting
def plot_repeats(repeat_class,bins,minval,maxval):
    global row_tracker; global column_tracker
    subdf=repeats.loc[repeats['repeat_class']==repeat_class]
    subdf=subdf[subdf['piq_begin']<5040]
    upstream_repeats=subdf['piq_begin']#upstream repeats
    downstream=downstream_repeats.loc[downstream_repeats['repeat_class']==repeat_class]['dpiq_begin']#downstream repeats
    full_rep=pd.concat([upstream_repeats,downstream])
    full_rep=full_rep[(full_rep>minval)&(full_rep<maxval)]
    kwargs={'edgecolor':'black','lw':0.05}
    axs[row_tracker,column_tracker].hist(full_rep,bins=bins,density=True,stacked=True,color='orange',**kwargs)
    axs[row_tracker,column_tracker].fill_between([4960,5040],axs[row_tracker,column_tracker].get_ylim()[0],axs[row_tracker,column_tracker].get_ylim()[1],alpha=.5,color='lightblue')
    axs[row_tracker,column_tracker].set_ylim(0,0.02)
    #title begins
    if '_' in repeat_class:
        axs[row_tracker,column_tracker].set_title(repeat_class.replace('_','\n'),fontsize=12.5)
    elif '/' in repeat_class:
        axs[row_tracker,column_tracker].set_title(repeat_class.replace('/','\n'),fontsize=12.5)
    else:
        axs[row_tracker,column_tracker].set_title(repeat_class,fontsize=15)
    #title ends
    axs[row_tracker,column_tracker].set_xticks([minval,5000,maxval])
    axs[row_tracker,column_tracker].set_xticklabels(['-'+str(5000-minval),'numt','+'+str(maxval-5000)],rotation=45,fontsize=12)
    axs[row_tracker,column_tracker].set_yticks([0,0.005,0.01,0.015,0.02])
    axs[row_tracker,column_tracker].set_yticklabels([0,'',0.01,'',0.02],fontsize=12)
    column_tracker+=1
    if column_tracker==4:
        column_tracker+=(-column_tracker)
        row_tracker+=1

#the strange repeats
sig_repeats=pd.Series(['Simple_repeat','Low_complexity', 'SINE/MIR','LINE/L2',
                       'LTR/ERVL', 'DNA/hAT-Charlie', 'LTR/ERV1','tRNA',
                       'LINE/L1','LTR/ERVL-MaLR','SINE/Alu','SINE/B4',
                      'DNA/hAT-Tip100','LINE/CR1','LTR/ERVK','SINE/B2'])

#visualize results and save the figure
row_tracker=0
column_tracker=0
fig,axs=plt.subplots(4,4,figsize=(8,6),sharex=True,sharey=True)
plt.subplots_adjust(wspace=0.1,
                    hspace=0.4)
sig_repeats.apply(plot_repeats,args=(100,4800,5200,))
plt.tight_layout(pad=1, w_pad=0.0001, h_pad=.1)
plt.text(-.55,-.55,'Density',fontsize=16,transform=ax[1,0].transAxes,rotation=90)
plt.savefig('../results/rm_densities.png',bbox_inches='tight',dpi=400)