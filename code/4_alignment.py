#import the required modules
import os
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

#function for creating LASTAL database for each organisms and align reversed mt dna with nuclear genome
def align_sequences(organism_name):
    os.chdir(f'../data/{organism_name}/')#change the working directory for first organism
    print(f'Building LASTAL database for {organism_name}...')
    call('lastdb db genome.fa', shell=True)#building database
    print(f'LASTAL database for {organism_name} has been successfully built!')
    call('lastal -r1 -q1 -a7 -b1 db d_mt.fa > d_mt_alignment.fa', shell=True)#align the genome and double mt dna into a file called d_mt_alignment.fa
    os.chdir('../../code/')#change back to the default 'code/' directory
    
organisms.apply(align_sequences)

#function for getting the e-value threshold and mask the significant alignments based on that value
def signifcant_alignments(organism_name):
    e_threshold=10**-4
    with open(f'../data/{organism_name}/d_mt_alignment.fa')as infile, open(f'../results/{organism_name}_significant_alignments.fa','w')as outfile:
        content=infile.readlines()
        print(f'Searching significant alignments within the {organism_name} genome...')
        for index, line in enumerate(content):
            if 'score' in line:
                e_value=float(line.rsplit()[3].split('=')[1])
                g_sequence = content[index + 1]
                mt_sequence = content[index + 2]
                if e_value < e_threshold:
                    outfile.write(line)
                    outfile.write(g_sequence)
                    outfile.write(mt_sequence)
                    outfile.write('\n')
        print(f'Significant alignments within the {organism_name} genome has been found!')
        print(f'Removing {organism_name} specific directory...')
        os.rmdir(f'../data/{organism_name}')#remove the organism specific directory not to have that much data
        print(f'{organism_name} specific data has been removed!')

organisms.apply(signifcant_alignments)