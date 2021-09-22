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

#function for getting genome sequence
def get_genome(organism_name):
    global problematic_genomes
    try:
        #get the latest genome version
        output_dir=os.path.join(f'../data/{organism_name}/')
        filename='CHECKSUMS'
        output_filepath=output_dir+filename
        checksums_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{filename}'
        urllib.request.urlretrieve(checksums_url, output_filepath)
        genome_version=''
        #get the genome version
        with open(f'../data/{organism_name}/CHECKSUMS')as infile:
            content=pd.Series(infile.readlines())
            mask=content.apply(lambda line: organism_name in line)
            genome_version=list(content[mask].apply(lambda line: line.split(organism_name+'.')[1].split('.dna')[0]))[0]
        #download genome sequence based on the previously defined genome version
        genome_name=f'{organism_name}.{genome_version}.dna.toplevel.fa.gz'
        genome_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{genome_name}'
        genome_filepath=output_dir+genome_name
        urllib.request.urlretrieve(genome_url, genome_filepath)
        #decompress gunzipped genome sequence
        with gzip.open(genome_filepath, 'rb')as infile, open(f'../data/{organism_name}/{genome_name[:-3]}', 'wb')as outfile:
            shutil.copyfileobj(infile, outfile)
        #remove the gunzipped genome
        os.remove(genome_filepath)
    except:
        problematic_genomes.append(organism_name)
    
#function for getting mt sequence
def get_mt(organism_name):
    #get the latest genome version
    output_dir=os.path.join(f'../data/{organism_name}/')
    filename='CHECKSUMS'
    output_filepath=output_dir+filename
    checksums_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{filename}'
    urllib.request.urlretrieve(checksums_url, output_filepath)
    genome_version=''
    #get the genome version
    with open(f'../data/{organism_name}/CHECKSUMS')as infile:
        content=pd.Series(infile.readlines())
        mask=content.apply(lambda line: organism_name in line)
        genome_version=list(content[mask].apply(lambda line: line.split(organism_name+'.')[1].split('.dna')[0]))[0]
    #download mt sequence based on the previously defined genome version
    mt_name=f'{organism_name}.{genome_version}.dna.chromosome.MT.fa.gz'
    mt_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{mt_name}'
    mt_filepath=output_dir+mt_name
    urllib.request.urlretrieve(mt_url, mt_filepath)
    #decompress gunzipped mt sequence
    with gzip.open(mt_filepath, 'rb')as infile, open(f'../data/{organism_name}/{mt_name[:-3]}', 'wb')as outfile:
        shutil.copyfileobj(infile, outfile) 
    #remove the gunzipped mt
    os.remove(mt_filepath)
    
