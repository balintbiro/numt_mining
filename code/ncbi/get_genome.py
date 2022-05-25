#import dependencies
import os
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
            call(f'wget --output-document=../../data/genomes/gDNA.fna.gz {URL+filename}',shell=True)
            #decompress genome
            call(f'gzip -d ../../data/genomes/gDNA.fna.gz',shell=True)
        except:
            print('A problem occured during gDNA acquisition!\nPossibilities are:\n\t-broken FTP connection\n\t-non existing assembly version')

    def get_mtDNA(self):
        try:
            #get mitochondrial id
            mtID=run(
                f'egrep \"{self.organism_name.replace('_',' ')} mitochondrion\" ../../data/genomes/mitochondrion.1.1.genomic.fna',
                shell=True,
                capture_output=True
            )
            mtID=str(mtID.stdout).split()[0][3:]
            #get mitochondrial sequence
            call(f'samtools faidx ../../data/genomes/mitochondrion.1.1.genomic.fna {mtID} > ../../data/genomes/mtDNA.fna',shell=True)
            #get duplicated mitochondria
            mtRecord=SeqIO.read("../../data/genomes/mtDNA.fna", "fasta")
            mtSeq=str(mtRecord.seq)
            with open('../../data/genomes/dmtDNA.fna','w')as outfile:
                outfile.write(str(mtID)+'\n'+2*mtSeq)
        except:
            print('A problem occured during mtDNA acquisition!\nPossibilities are:\n\t-no available mtDNA file for the given organism')

    def LASTAL(self):
        try:
            #generate LASTAL db
            call('lastdb ../../data/genomes/db ../../data/genomes/gDNA.fna', shell=True)
            #align gDNA with dmtDNA
            call('lastal -r1 -q1 -a7 -b1 ../../data/genomes/db dmtDNA.fna -f BlastTab > dmtDNA.tab', shell=True)
        except:
            print('A problem occured during db building or alignment!\nPossibilities are:-gDNA and/or mtDNA is/are missing')


#get mitochondrion file 
call(
    f'wget --directory-prefix=../../data/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz',
    shell=True
)

#decompress mitochondrion file
call(
    f'gzip -d ../../data/genomes/mitochondrion.1.1.genomic.fna.gz',
    shell=True
)

#connect to NCBI FTP site
ftp=FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()
#change directory to mammalian genomes
ftp.cwd('/genomes/refseq/vertebrate_mammalian/')
#get organism names
organisms=ftp.nlst()[:1]

#create folder for genomes
if os.path.isdir('../../data/genomes/')==False:
    os.mkdir('../../data/genomes/')

#create classes and apply functions
for organism in organisms:
    organism_class=process_DNA(organism)#instantiation
    organism_class.get_gDNA()#download and decompress gDNA sequence
    organism_class.get_mtDNA()#get mtDNA sequence
    organism_class.LASTAL()#building LASTAL db and align dmtDNA and gDNA