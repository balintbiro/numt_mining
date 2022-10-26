#import dependencies
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from ftplib import FTP
from subprocess import call, run

#class definition for processing DNA
class process_DNA():
    def __init__(self,organism_name,output_file):
        self.organism_name=organism_name.capitalize()
        self.output_file=output_file

    def get_gDNA(self):
        try:
            if os.path.getsize(f'../data/alignments/{self.organism_name}_numts.csv')>1000:
                #connect to the NCBI FTP site
                ftp=FTP('ftp.ncbi.nlm.nih.gov')
                ftp.login()
                #go to the latest assembly version of the given organism
                ftp.cwd(f'/genomes/refseq/vertebrate_mammalian/{self.organism_name}/latest_assembly_versions/')
                LatestAssembly=ftp.nlst()[0]
                URL=f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/{self.organism_name}/latest_assembly_versions/{LatestAssembly}/'
                filename=f'{LatestAssembly}_genomic.fna.gz'
                #download latest genome version
                call(f'wget --output-document=../data/genomes/gDNA.fna.gz {URL+filename}',shell=True)
                #decompress genome
                call(f'gzip -d ../data/genomes/gDNA.fna.gz',shell=True)
        except:
            print(f'A problem occured during {self.organism_name} gDNA acquisition!')

    def get_random_seq(self):
        try:
            alignments=pd.read_csv(f'../data/alignments/{self.organism_name}_numts.csv')
            numt_number=len(alignments)
            #get min and max lengths of numts
            numt_range=alignments['genomic_length']
            record_dict=SeqIO.to_dict(SeqIO.parse("../data/genomes/gDNA.fna", "fasta"))
            genomic_ids=pd.Series(record_dict.keys())
            #random sampling of genomic ids
            np.random.seed(0)
            genomic_ids=pd.Series(np.random.choice(genomic_ids,size=int(numt_number*1.5),replace=True))
            #get the corresponding sequences
            sequences=genomic_ids.apply(lambda g_id: str(record_dict[g_id].seq))
            #get the corresponding sequence lengths
            sizes=genomic_ids.apply(lambda g_id: len(str(record_dict[g_id].seq)))
            #create df from selected ids
            df=pd.DataFrame(data=[genomic_ids,sizes,sequences]).T
            df.columns=['genomic_id','size','sequence']
            #add starting positions
            df['sample_start']=df.apply(lambda row: np.random.choice(row['size']),axis=1)
            #create inner function for sampling
            def random_sample(row):
                start=row['sample_start']-5000
                np.random.seed(int(row['size']))
                sample_size=np.random.choice(numt_range)
                end=row['sample_start']+sample_size+5000
                if start<0:
                    start=0
                if end>len(row['sequence']):
                    end=len(row['sequence'])
                sample_seq=row['sequence'][start:end]
                upstream_size=row['sample_start']-start
                downstream_size=end-(row['sample_start']+sample_size)
                random_sequences.append([self.organism_name,row['genomic_id'],row['size'],row['sample_start'],upstream_size,sample_seq,sample_size,downstream_size,'random'])
            #global variable
            random_sequences=[]
            #make the function work
            df.apply(random_sample,axis=1)
            random_sequences=pd.DataFrame(
                data=random_sequences,
                columns=['organism_name','genomic_id','genomic_size','genomic_start','upstream_size','sequence','sample_size','downstream_size','label']
            )
            random_sequences.to_csv(self.output_file,mode='a',index=False,header=False)
        except:
            print(f'Problem with {self.organism_name} during random sampling!')


    def clear_folder(self):
        call('rm ../data/genomes/*',shell=True)

#create folder for genomes
if os.path.isdir('../data/genomes/')==False:
    os.mkdir('../data/genomes/')

#create empty output file and write it
ml_df=pd.DataFrame(
		data=[],
		columns=['organism_name','genomic_id','genomic_size','genomic_start','upstream_size','sequence','sample_size','downstream_size','label']
	)
ml_df.to_csv('../data/ml_input_randoms.csv',index=False)

#read in numts
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#get organism_names
organisms=numts['organism_name'].drop_duplicates().tolist()

#create classes and apply functions
for organism in organisms:
    if os.path.getsize('../data/ml_input_randoms.csv')<10000000000:#size is under 10 GB just to make sure we dont have an overflow
        organism_class=process_DNA(organism,'../data/ml_input_randoms.csv')#instantiation
        organism_class.get_gDNA()#download and decompress gDNA sequence
        organism_class.get_random_seq()#create random sequence sample
        organism_class.clear_folder()#delete the unnecessary files