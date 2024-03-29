#importing the required modules
import os
import pandas as pd
from subprocess import call

#downloading the genome summary file
if os.path.exists('../data/assembly_summary_genbank.txt')==False:
    call('wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -P ../data/', shell=True)

#add header from assembly summary file to genomic urls files
call("awk 'NR==2' ../data/assembly_summary_genbank.txt > ../data/genomic_urls.tsv", shell=True)

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
organisms.apply(lambda organism_name: call(f"grep -E '{organism_name}' ../data/assembly_summary_genbank.txt | grep 'reference genome\|representative genome'>> ../data/genomic_urls.tsv", shell=True))

#load the tsv files which contains the selected genomes into a pandas dataframe
df=pd.read_csv('../data/genomic_urls.tsv',sep='\t',header=0)
df=df.set_index('organism_name')

#get the genome ftp links
ftp_paths=df['ftp_path']

#defining a function for downloading genomes
def get_genome(ftp_path):
    organism_name=ftp_paths[ftp_paths==ftp_path].index[0].replace(' ','_')
    filename=ftp_path.split('/')[-1]+'_genomic.fna.gz'
    if os.path.exists(f'../data/{organism_name}/')==False:
        os.mkdir(f'../data/{organism_name}/')
    url=ftp_path+'/'+filename
    call(f'wget {url} -P ../data/{organism_name}/', shell=True)
    call(f'gunzip -c ../data/{organism_name}/{filename} > ../data/{organism_name}/genome.fa', shell=True)
    
ftp_paths.apply(get_genome)