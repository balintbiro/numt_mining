#importing the required modules
import os
import pandas as pd

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(','))
organisms.apply(mkdir)

#get a merged list for all the informations about alignments
features = []
with open ('../data/significant_alignments.fasta') as infile:
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
            features.append([score, eg2_value, e_value, genomic_id, genomic_start,
                            mitochondrial_start, genomic_length, mitochondrial_length, genomic_strand,
                            mitochondrial_strand, genomic_size, genomic_sequence,
                             mitochondrial_sequence])
                         
