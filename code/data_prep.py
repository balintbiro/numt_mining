#importing the required modules
import os
import gzip
import shutil
import pandas as pd
import urllib.request   

#function for creating repository for every organisms
def mkdir(organism_name):
    parent_dir='../data/'
    path=os.path.join(parent_dir+organism_name)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    
organisms=pd.Series(['Oryctolagus_cuniculus'])
organisms.apply(mkdir)

#function for getting the latest genome version and required sequences
def get_data(organism_name):
    #get the latest genome version
    output_dir=os.path.join(f'../data/{organism_name}/')
    filename='CHECKSUMS'
    output_filepath=output_dir+filename
    checksums_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{filename}'
    urllib.request.urlretrieve(checksums_url, output_filepath)
    genome_version=''
    #get the sequences based on the previously acquired genome version
    with open(f'../data/{organism_name}/CHECKSUMS')as infile:
        content=pd.Series(infile.readlines())
        mask=content.apply(lambda line: organism_name in line)
        genome_version=list(content[mask].apply(lambda line: line.split(organism_name+'.')[1].split('.dna')[0]))[0]
    genome_name=f'{organism_name}.{genome_version}.dna.toplevel.fa.gz'
    genome_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{genome_name}'
    genome_filepath=output_dir+genome_name
    urllib.request.urlretrieve(genome_url, genome_filepath)
    with gzip.open(genome_filepath, 'rb')as infile, open(f'../data/{organism_name}/{genome_name[:-3]}', 'wb')as outfile:
        shutil.copyfileobj(infile, outfile)
    
organisms.apply(get_data)