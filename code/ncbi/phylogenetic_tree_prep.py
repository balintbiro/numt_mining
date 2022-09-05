#import dependencies
import os
import pandas as pd
from subprocess import call, run

#read numts df
numts=pd.read_csv('../data/ncbi_numts_p26.csv')

#get organism_names
organism_names=pd.Series(numts['organism_name'].unique()).apply(lambda name: name[0].upper()+name[1:])

#function for extracting mtDNA sequences
def get_mtDNA(organism_name):
    global organisms
    global mtIDs
    try:
        #get mitochondrial id
        mtID=run(f"""egrep '{organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.1.1.genomic.fna | grep mitochondrion""",shell=True,capture_output=True)
        mtID=str(mtID.stdout).split()[0][3:]
        #create file for mtDNAs
        call('touch ../data/genomes/mtDNAs_for_phylogenetics.fna',shell=True)
        #get previous mtDNAs file size
        prev_size=os.path.getsize('../data/genomes/mtDNAs_for_phylogenetics.fna')
        #get mitochondrial sequence
        call(f'samtools faidx ../data/genomes/mitochondrion.1.1.genomic.fna {mtID} >> ../data/genomes/mtDNAs_for_phylogenetics.fna',shell=True)
        #some organisms are not present in the mitochondrion.1.1.genomic.fna file so try mitochondrion.2.1.genomic.fna instead
        if (os.path.getsize('../data/genomes/mtDNAs_for_phylogenetics.fna')-prev_size)<1000:
            mtID=run(f"""egrep '{organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.2.1.genomic.fna | grep mitochondrion""",shell=True,capture_output=True)
            mtID=str(mtID.stdout).split()[0][3:]
            call(f'samtools faidx ../data/genomes/mitochondrion.2.1.genomic.fna {mtID} >> ../data/genomes/mtDNAs_for_phylogenetics.fna',shell=True)
        organisms.append(organism_name)
        mtIDs.append(mtID)
    except:
        print(f'No sequence for {organism_name}!')

#apply the function
organisms=[]
mtIDs=[]
organism_names.apply(get_mtDNA)

#create pd series from mtids and organism names
mtDNAs=pd.Series(
        data=organisms,
        index=mtIDs
    )

#write series into csv file
mtDNAs.to_csv('../data/mtDNA_ids.csv')

#create multiple sequence alignment
call('clustalo --infile=../data/genomes/mtDNAs_for_phylogenetics.fna --seqtype=DNA --outfile=../results/aligned_mtDNAs.phy',shell=True)