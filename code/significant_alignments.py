#import the required modules
import os
import pandas as pd

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(','))
organisms.apply(mkdir)

#function for getting the e-value threshold and mask the significant alignments based on that value
def signifcant_alignments(organism_name):
    home_dir=os.path.join(f'../data/{organism_name}/')
    e_values=[]
    with open(home_dir+'r_mt_alignment.fa')as infile:
        content=pd.Series(infile.readlines())
        mask=content.apply(lambda line: 'EG2' in line)
        content[mask].apply(lambda line: e_values.append(float(line.rsplit()[3].split('=')[1])))
    e_values.sort()
    e_threshold=e_values[0]
    with open(home_dir+'d_mt_alignment.fa')as infile, open(f'../results/{organism_name}_signifcant_alignments.fa','w')as outfile:
        content=infile.readlines()
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
                
organisms.apply(signifcant_alignments)