#import the required modules
import os
import pandas as pd
from subprocess import call

#function for writing reversed mt sequence
def create_r_mt(organism_name):
    organism_dir=os.path.join(f'../data/{organism_name}/')
    mt_filename=f'{organism_name}.mt.fa'
    mt_path=organism_dir+mt_filename
    header=''
    mt_seq=''
    with open(mt_path)as infile, open(organism_dir+'r_mt_seq.fa','w')as outfile:
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