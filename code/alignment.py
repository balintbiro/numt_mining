#import the required modules
import os
import pandas as pd
from subprocess import call

#function for writing reversed mt sequence
def create_r_mt(organism_name):
    organism_dir=os.path.join(f'../data/{organism_name}/')
    mt_filename='mt.fa'
    mt_path=organism_dir+mt_filename
    header=''
    mt_seq=''
    with open(mt_path)as infile, open(organism_dir+'r_mt.fa','w')as outfile:
        for line in infile:
            if '>' in line:
                header=line
            else:
                mt_seq+=line[:-1]
        r_mt_seq=mt_seq[::-1]
        outfile.write(header)
        outfile.write(r_mt_seq)

organisms=pd.Series(['Oryctolagus_cuniculus'])


organisms.apply(create_r_mt)

#function for creating LASTAL database for each organisms and align reversed mt dna with nuclear genome
def align_r_mt(organism_name):
    os.chdir(f'../data/{organism_name}')#change the working directory
    call('lastdb db genome.fa', shell=True)#building database
    call('lastal db r_mt.fa > r_mt_alignment.fa', shell=True)#align the genome and reversed mt dna into a file called r_mt_alignment.fa
    
organisms.apply(align_r_mt)