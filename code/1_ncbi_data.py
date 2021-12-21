#importing the required modules
import os
import pandas as pd
from subprocess import call

#downloading the genome summary file
if os.path.exists('../data/assembly_summary_genbank.txt')==False:
    call('wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -P ../data/', shell=True)

#function for creating repository for every organisms
def mkdir(organism_name):
    parent_dir='../data/'
    path=os.path.join(parent_dir+organism_name)
    if os.path.exists(path)==False:
        os.mkdir(path)

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(',')).apply(lambda organism_name: organism_name.replace('_',' '))

#get the downloadable urls for the required organisms
print('Get the downloadable URLs of the genomes for the required organisms...')
organisms.apply(lambda organism_name: call(f"grep -E '{organism_name}' ../data/assembly_summary_genbank.txt >> ../data/genomic_urls.tsv", shell=True))