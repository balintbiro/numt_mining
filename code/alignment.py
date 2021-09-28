#import the required modules
import os
import pandas as pd
from subprocess import call

#get the organism names
with open('../data/organism_names.txt')as infile:
    organisms=pd.Series(infile.readlines()[0].split(','))

#function for writing reversed and double mt sequence
def mt_versions(organism_name):
    organism_dir=os.path.join(f'../data/{organism_name}/')
    mt_filename='mt.fa'
    mt_path=organism_dir+mt_filename
    header=''
    mt_seq=''
    with open(mt_path)as infile, open(organism_dir+'r_mt.fa','w')as r_outfile, open(organism_dir+'d_mt.fa','w')as d_outfile:
        for line in infile:
            if '>' in line:
                header=line
            else:
                mt_seq+=line[:-1]
        r_mt_seq=mt_seq[::-1]
        r_outfile.write(header)
        r_outfile.write(r_mt_seq)
        d_mt_seq=2*mt_seq
        d_outfile.write(header)
        d_outfile.write(d_mt_seq)


organisms.apply(mt_versions)

#function for creating LASTAL database for each organisms and align reversed mt dna with nuclear genome
def align_sequences(organism_name):
    os.chdir(f'/home/birobalint/numt_mining/data/{organism_name}')#change the working directory
    call('lastdb db genome.fa', shell=True)#building database
    call('lastal db r_mt.fa > r_mt_alignment.fa', shell=True)#align the genome and reversed mt dna into a file called r_mt_alignment.fa
    call('lastal -r1 -q1 -a7 -b1 db d_mt.fa > d_mt_alignment.fa', shell=True)#align the genome and double mt dna into a file called d_mt_alignment.fa
    
organisms.apply(align_sequences)

#function for getting the e-value threshold and mask the significant alignments based on that value
def signifcant_alignments(organism_name):
    os.chdir(f'/home/birobalint/numt_mining/data/{organism_name}')#change the working directory
    e_values=[]
    with open(f'/home/birobalint/numt_mining/data/{organism_name}/r_mt_alignment.fa')as infile:
        content=pd.Series(infile.readlines())
        mask=content.apply(lambda line: 'EG2' in line)
        content[mask].apply(lambda line: e_values.append(float(line.rsplit()[3].split('=')[1])))
    e_values.sort()
    e_threshold=e_values[0]
    with open(f'/home/birobalint/numt_mining/data/{organism_name}/d_mt_alignment.fa')as infile, open(f'/home/birobalint/numt_mining/results/{organism_name}_signifcant_alignments.fa','w')as outfile:
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