#import the required modules
import os
import pandas as pd
from subprocess import call

#function for writing reversed and double mt sequence
def mt_versions(organism_name):
    organism_dir=os.path.join(f'../data/{organism_name}/')
    mt_filename='mt.fa'
    mt_path=organism_dir+mt_filename
    header=''
    mt_seq=''
    with open(mt_path)as infile, open(organism_dir+'r_mt.fa','w')as r_outfile, open(organism_dir+'d_mt.fa','w')as d_outfile:
        for line in infile:
            if '>' in line:
                header=line
            else:
                mt_seq+=line[:-1]
        r_mt_seq=mt_seq[::-1]
        r_outfile.write(header)
        r_outfile.write(r_mt_seq)
        d_mt_seq=2*mt_seq
        d_outfile.write(header)
        d_outfile.write(d_mt_seq)

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(','))


organisms.apply(mt_versions)

#function for creating LASTAL database for each organisms and align reversed mt dna with nuclear genome
def align_sequences(organism_name):
    os.chdir(f'../data/{organism_name}')#change the working directory
    call('lastdb db genome.fa', shell=True)#building database
    call('lastal db r_mt.fa > r_mt_alignment.fa', shell=True)#align the genome and reversed mt dna into a file called r_mt_alignment.fa
    call('lastal -r1 -q1 -a7 -b1 db d_mt.fa > d_mt_alignment.fa', shell=True)#align the genome and double mt dna into a file called d_mt_alignment.fa
    
organisms.apply(align_sequences)