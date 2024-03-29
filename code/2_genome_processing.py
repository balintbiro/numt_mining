#import the required modules
import os
import filecmp
import pandas as pd
from subprocess import call

#get the original organism names
with open('../data/organism_names.txt')as infile:
    default_organisms=pd.Series(infile.readlines()[0].split(','))

#get the organism names that have problem(s)
problems_filepath='../data/problematic_organisms.txt'
problematic_organisms=[]
if os.path.exists(problems_filepath):
    with open(problems_filepath)as infile:
        problematic_organisms=pd.Series(infile.readlines()[0].split(','))
    
#get the organisms that have mt and g genomes too
organisms=pd.Series(list(set(default_organisms)-set(problematic_organisms)))

#define function for genome processing
def process_genome(organism_name):
    try:
        os.chdir(f'../data/{organism_name}/')#change the working directory for the first organism
    except FileNotFoundError:
        os.chdir(f'../{organism_name}/')#change the working directory for every other organisms
    call("""grep -E '>' genome.fa > headers.txt""", shell=True)#get all the headers from the genome
    call("""grep -E '>*REF' headers.txt > ref_headers.txt""", shell=True)#get the reference headers
    if filecmp.cmp('headers.txt','ref_headers.txt')==True:
        print(f'There are no alternative contigs in the nuclear genome of {organism_name}!')
    else:#that means that there are alternative headers
        print(f'The nuclear genome of {organism_name} contains alternative contigs! These contigs need to be removed for the alignment. This will take some time!')
        os.rename('genome.fa','alt_genome.fa')#rename the original genome.fa file to alt_genome.fa which contains the alternative sequences
        
        with open('alt_genome.fa')as infile, open('processed_genome.fa','w')as outfile: #transform alt_genome to one liner form
            for line in infile:
                if '>' in line:
                    outfile.write('\n')
                    outfile.write(line)
                else:
                    outfile.write(line[:-1])
        
        call("""grep -f ref_headers.txt -A1 processed_genome.fa > genome.fa""", shell=True)#get the reference sequences based on the reference headers in the ref_headers.txt file
        call("""sed -i '/^-/d' genome.fa""", shell=True)##delete the line starting with '-' sign which is is introduced by grep command
        call("""rm alt_genome.fa processed_genome.fa""", shell=True)#remove the unnecessary genome files
        print(f'The alternative contigs from the genome of {organism_name} were successfully removed! The files are ready for LASTAL.')
    os.chdir('../../code/')#change back to the default 'code/' directory
    
organisms.apply(process_genome)