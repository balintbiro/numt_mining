#importing the required modules
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

#define function for creating a pandas dataframe from the signifcant alignments
def create_df(organism_name):
    #get a merged list for all the informations about alignments
    features = []
    with open (f'../results/{organism_name}_significant_alignments.fa') as infile:
        content = infile.readlines()
        for index, line in enumerate(content):
            if 'score' in line:
                score = int(line.rsplit()[1].split('=')[1])
                eg2_value = float(line.rsplit()[2].split('=')[1])
                e_value = float(line.rsplit()[3].split('=')[1])
                
                #genomic informations
                genomic = content[index + 1]
                genomic_id = genomic.rsplit()[1]
                genomic_start = int(genomic.rsplit()[2])
                genomic_length = int(genomic.rsplit()[3])
                genomic_strand = genomic.rsplit()[4]
                genomic_size = int(genomic.rsplit()[5])
                genomic_sequence = genomic.rsplit()[6]
                
                #mitochondrial informations
                mitochondrial = content[index + 2]
                mitochondrial_start = int(mitochondrial.rsplit()[2])
                mitochondrial_length = int(mitochondrial.rsplit()[3])
                mitochondrial_strand = mitochondrial.rsplit()[4]
                mitochondrial_sequence = mitochondrial.rsplit()[6]
                
                #merge informations together
                features.append([score, eg2_value, e_value, genomic_id, genomic_start,
                                mitochondrial_start, genomic_length, mitochondrial_length, genomic_strand,
                                mitochondrial_strand, genomic_size, genomic_sequence,
                                 mitochondrial_sequence])
    #create dataframe and name columns
    alignment_df = pd.DataFrame(features)
    alignment_df.columns = ['score', 'eg2_value', 'e_value', 'g_id', 'g_start', 'mt_start',
                            'g_length', 'mt_length', 'g_strand', 'mt_strand', 'g_size',
                            'g_sequence', 'mt_sequence']
    #drop the duplicates where the genomic start ('g_start') and genomic length ('g_length') are equal but the mitochondrial start differs with the mitochondrial length
    #that means that we have the alignments from the duplicated mitochondria
    updated_alignments = alignment_df.drop_duplicates(subset = ['g_start', 'g_length'])
    
    #filter out the line where the alignment found the mt itself
    mt_mask=updated_alignments['g_length']!=updated_alignments['g_size']
    updated_alignments=updated_alignments[mt_mask]
    
    #write the dataframe into a csv file
    updated_alignments.to_csv(f'../results/{organism_name}_numts.csv', index = False)
    
organisms.apply(create_df)