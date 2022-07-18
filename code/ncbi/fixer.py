#import dependencies
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from ftplib import FTP
from subprocess import call, run

#class definition for processing DNA
class process_DNA():
    def __init__(self, organism_name):
        self.organism_name=organism_name

    def get_gDNA(self):
        try:
            ftp=FTP('ftp.ncbi.nlm.nih.gov')
            ftp.login()
            ftp.cwd(f'/genomes/refseq/vertebrate_mammalian/{self.organism_name}/latest_assembly_versions/')
            LatestAssembly=ftp.nlst()[0]
            URL=f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/{self.organism_name}/latest_assembly_versions/{LatestAssembly}/'
            filename=f'{LatestAssembly}_genomic.fna.gz'
            #download genome
            call(f'wget --output-document=../data/genomes/gDNA.fna.gz {URL+filename}',shell=True)
            #decompress genome
            call(f'gzip -d ../data/genomes/gDNA.fna.gz',shell=True)
        except:
            print(f'A problem occured during {self.organism_name} gDNA acquisition!\nPossibilities are:\n\t-broken FTP connection\n\t-non existing assembly version')

    def get_mtDNA(self):
        try:
            #get mitochondrial id
            mtID=run(f"""egrep '{self.organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.1.1.genomic.fna | grep mitochondrion""",shell=True,capture_output=True)
            mtID=str(mtID.stdout).split()[0][3:]
            #get mitochondrial sequence
            call(f'samtools faidx ../data/genomes/mitochondrion.1.1.genomic.fna {mtID} > ../data/genomes/mtDNA.fna',shell=True)
            if os.path.getsize('../data/genomes/mtDNA.fna')<1000:
                mtID=run(f"""egrep '{self.organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.1.1.genomic.fna | grep mitochondrion""",shell=True,capture_output=True)
                mtID=str(mtID.stdout).split()[0][3:]
                call(f'samtools faidx ../data/genomes/mitochondrion.2.1.genomic.fna {mtID} > ../data/genomes/mtDNA.fna',shell=True)
            #get duplicated mitochondria
            mtRecord=SeqIO.read("../data/genomes/mtDNA.fna", "fasta")
            mtSeq=str(mtRecord.seq)
            with open('../data/genomes/dmtDNA.fna','w')as outfile:
                outfile.write('>'+str(mtID)+'\n'+2*mtSeq)
        except:
            print(f'A problem occured during {self.organism_name} mtDNA acquisition!\nPossibilities are:\n\t-no available mtDNA sequence for the given organism\n\t-problematic ID')

    def LASTAL(self):
        try:
            if os.path.getsize('../data/genomes/mtDNA.fna')>1000:
                #generate LASTAL db
                call('lastdb ../data/genomes/db ../data/genomes/gDNA.fna', shell=True)
                #align gDNA with dmtDNA
                call('lastal -r1 -q1 -a7 -b1 ../data/genomes/db ../data/genomes/dmtDNA.fna  > ../data/genomes/aligned_dmtDNA.afa', shell=True)
            else:
                print(f'A problem occured during {self.organism_name} db building or alignment!\nPossibilities are:\n\t-gDNA and/or mtDNA is/are missing')
        except:
            print(f'A problem occured during {self.organism_name} db building or alignment!\nPossibilities are:\n\t-gDNA and/or mtDNA is/are missing')

    def filter_alignments(self):
        try:
            #get data from alignment file
            lines = []
            with open ('../data/genomes/aligned_dmtDNA.afa') as infile:
                content = infile.readlines()
                for index, line in enumerate(content):
                    if 'score' in line:
                        score = int(line.rsplit()[1].split('=')[1])
                        eg2_value = float(line.rsplit()[2].split('=')[1])
                        e_value = float(line.rsplit()[3].split('=')[1])
                        #genomic informations
                        genomic = content[index + 1]
                        genomic_id,genomic_start,genomic_length,genomic_strand,genomic_size,genomic_sequence = genomic.rsplit()[1],int(genomic.rsplit()[2]),int(genomic.rsplit()[3]),genomic.rsplit()[4],int(genomic.rsplit()[5]),genomic.rsplit()[6] 
                        #mitochondrial informations
                        mitochondrial = content[index + 2]
                        mitochondrial_start,mitochondrial_length,mitochondrial_strand,mitochondrial_sequence = int(mitochondrial.rsplit()[2]),int(mitochondrial.rsplit()[3]),mitochondrial.rsplit()[4],mitochondrial.rsplit()[6]
                        lines.append([score, eg2_value, e_value, genomic_id, genomic_start,
                                        mitochondrial_start, genomic_length, mitochondrial_length, genomic_strand,
                                        mitochondrial_strand, genomic_size, genomic_sequence,
                                         mitochondrial_sequence])
            #create the alignment df
            alignments=pd.DataFrame(
                data=lines,
                columns=['score', 'eg2_value', 'e_value', 'genomic_id', 'genomic_start',
                        'mitochondrial_start', 'genomic_length', 'mitochondrial_length', 'genomic_strand',
                        'mitochondrial_strand', 'genomic_size', 'genomic_sequence',
                         'mitochondrial_sequence']
            )
            #create filter to discard artifacts that are the results of using dmtDNA
            mtRecord=SeqIO.read("../data/genomes/mtDNA.fna", "fasta")
            mtSize=len(str(mtRecord.seq))
            size_fil=alignments['mitochondrial_start']<mtSize
            #create filter to discard non significant alignments
            sig_fil=alignments['e_value']<10**-3
            #create filters to discard the mitochondrion itself (some toplevel fasta files contains also the mitochondrion)
            mt_fil=alignments['genomic_size'].apply(
                lambda genomic_size: genomic_size not in np.arange(mtSize-5,mtSize+5,1)
                )
            #apply filters and write output into alignments folder
            alignments=alignments[size_fil][sig_fil][mt_fil]
            alignments.to_csv(f'../data/alignments/{self.organism_name}_numts.csv',header=True)
        except:
            print(f'A problem occured during {self.organism_name} alignment csv construction!\nPossibilities are:\n\t-missing input files\n\t-no alignment')

    def clear_folder(self):
        folder_content=pd.Series(os.listdir('../data/genomes/'))
        necessary_files=['mitochondrion.1.1.genomic.fna','mitochondrion.1.1.genomic.fna.fai','mitochondrion.2.1.genomic.fna','mitochondrion.2.1.genomic.fna.fai']
        folder_content.apply(
            lambda file: os.remove(f'../data/genomes/{file}') if file not in necessary_files else file
        )

#get mitochondrion files and decompress them
#call(f'wget --directory-prefix=../data/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz',shell=True)
#call(f'wget --directory-prefix=../data/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.1.genomic.fna.gz',shell=True)
#call(f'gzip -d ../data/genomes/mitochondrion.1.1.genomic.fna.gz',shell=True)
#call(f'gzip -d ../data/genomes/mitochondrion.2.1.genomic.fna.gz',shell=True)

#connect to NCBI FTP site
ftp=FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
#change directory to mammalian genomes
ftp.cwd('/genomes/refseq/vertebrate_mammalian/')
#get organism names
organisms=ftp.nlst()

#already aligned
already_aligned=pd.Series(os.listdir('../data/alignments/')).apply(lambda filename: filename[:-10]).tolist()
organisms=list(set(organisms)-set(already_aligned))

#create folder for genomes
if os.path.isdir('../data/genomes/')==False:
    os.mkdir('../data/genomes/')

#create folder for significant alignments
if os.path.isdir('../data/alignments/')==False:
    os.mkdir('../data/alignments/')

#create classes and apply functions
for organism in organisms:
    organism_class=process_DNA(organism)#instantiation
    organism_class.get_gDNA()#download and decompress gDNA sequence
    organism_class.get_mtDNA()#get mtDNA sequence
    organism_class.LASTAL()#building LASTAL db and align dmtDNA and gDNA
    organism_class.filter_alignments()#filter alignments based on mt start, genomic size, e value
    organism_class.clear_folder()#delete the unnecessary files
