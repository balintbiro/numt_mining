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

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(','))
organisms.apply(mkdir)

#function for getting genome sequence
def get_genome(organism_name, genome):#genome can be nuclear or mitochondrial
    global problems
    #get the latest genome version
    output_dir=f'../data/{organism_name}/'
    filename='CHECKSUMS'
    output_filepath=os.path.join(output_dir+filename)
    checksums_url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{filename}'
    urllib.request.urlretrieve(checksums_url, output_filepath)
    genome_version=''
    #get the latest genome version
    with open(f'../data/{organism_name}/CHECKSUMS')as infile:
        content=pd.Series(infile.readlines())
    mt_mask=content.apply(lambda row: '.' in row
               and ('MT' in row.split('.')[-3] or 'mitochondr' in row.split('.')[-3])
              and 'dna' in row.split('.'))
    g_mask=content.apply(lambda row: '.' in row
              and 'dna' in row.split('.')
                and 'toplevel' in row.split('.'))
    #download genome sequence based on the previously defined genome version
    if genome=='nuclear':
        filename=content[g_mask].tolist()[0].rsplit()[-1]
        filename_short='genome.fa'
    elif genome=='mitochondrial':
        filename=content[mt_mask].tolist()[0].rsplit()[-1]
        filename_short='mt.fa'
    try:
        url=f'http://ftp.ensembl.org/pub/release-104/fasta/{organism_name.lower()}/dna/{filename}'
        filepath=output_dir+filename
        urllib.request.urlretrieve(url, filepath)
        #decompress gunzipped genome sequence
        with gzip.open(filepath, 'rb')as infile, open(f'../data/{organism_name}/{filename_short}', 'wb')as outfile:
            shutil.copyfileobj(infile, outfile)
        #remove the gunzipped genome
        os.remove(filepath)
    except:
        problems.append(organism_name)

#function for writing problematic organisms
def handling_exceptions(problem, genome):
    if len(problems)==0:
        pass
    else:
        report_dir='../results/problem_reports/'
        filename=''
        if genome=='nuclear':
            filename='problematic_genomes.txt'
        elif genome=='mitochondrial':
            filename='problematic_mts.txt'
        if os.path.exists(report_dir):
            with open(os.path.join(report_dir)+filename,'w')as output:
                output.write(f'No available {genome} genome for {problem}\n')
        else:
            os.mkdir(report_dir)
            with open(os.path.join(report_dir)+filename,'w')as output:
                output.write(f'No available {genome} genome for {problem}\n')

#global variable for nuclear problems      
problems=[]
organisms.apply(get_genome,args=('nuclear',))
problems=pd.Series(problems)
problems.apply(handling_exceptions, args=('nuclear',))

#global variable for mitochondrial problems      
problems=[]
organisms.apply(get_genome,args=('mitochondrial',))
problems=pd.Series(problems)
problems.apply(handling_exceptions, args=('mitochondrial',))