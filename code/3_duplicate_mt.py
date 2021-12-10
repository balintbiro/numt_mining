#import the required modules
import os
import pandas as pd

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

#function for duplicate mt sequence
def duplicate_mt(organism_name):
    organism_dir=os.path.join(f'../data/{organism_name}/')
    mt_filename='mt.fa'
    mt_path=organism_dir+mt_filename
    with open(mt_path)as infile, open(organism_dir+'d_mt.fa','w')as d_outfile:
        content=pd.Series(infile.readlines())
        header_mask=content.apply(lambda line: '>' in line)
        sequence_mask=content.apply(lambda line: '>' not in line)
        sequence_parts=[]
        content[sequence_mask].apply(lambda line: sequence_parts.append(line[:-1]))
        sequence=''.join(sequence_parts)
        duplicated_sequence=2*sequence
        header=content[header_mask].tolist()[0]
        d_outfile.write(header)
        d_outfile.write(duplicated_sequence)


organisms.apply(duplicate_mt)