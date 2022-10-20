#import dependencies
import pandas as pd

#read in dataframes
numts=pd.read_csv('../data/ncbi_numts_p26.csv')
repeatranks=pd.read_csv('../results/repeatranks.csv')

#add repeatranks 1st to numts
numts['u_1st_repeatclass']=repeatranks['u_1st_repeatclass']
numts['d_1st_repeatclass']=repeatranks['d_1st_repeatclass']

#create filters and filter numts df
ufil=numts['upstream_5kb'].apply(lambda seq:type(seq)==str and len(seq)>10)
dfil=numts['downstream_5kb'].apply(lambda seq:type(seq)==str and len(seq)>10)
numts=numts[ufil][dfil]
numts=numts[['upstream_5kb','genomic_sequence','downstream_5kb','order','u_1st_repeatclass','d_1st_repeatclass']]

#get seqlogo input aka 10 bp upstream flanking, 10 bp numt upstream end, 10 bp numt downstream end, 10 bp downstream flanking
numts['seqlogo_input']=numts.apply(lambda row: row['upstream_5kb'][-10:]+row['genomic_sequence'][:10]+row['genomic_sequence'][-10:]+row['downstream_5kb'][:10],axis=1)
numts=numts[['seqlogo_input','order','u_1st_repeatclass','d_1st_repeatclass']]

#write output
numts.to_csv('../data/seqlogo_inputs.csv')