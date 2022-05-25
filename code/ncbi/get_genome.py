#import dependencies
import os
import pandas as pd
from Bio import SeqIO
from ftplib import FTP
from subprocess import call, run

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
organisms=pd.Series(ftp.nlst()[:1])

#create folder for genomes
if os.path.isdir('../../data/genomes/')==False:
    os.mkdir('../../data/genomes/')

#function for dealing with the corresponding genomes
def process_DNA(organism_name):
    ftp=FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd(f'/genomes/refseq/vertebrate_mammalian/{organism_name}/latest_assembly_versions/')
    LatestAssembly=ftp.nlst()[0]
    URL=f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/{organism_name}/latest_assembly_versions/{LatestAssembly}/'
    filename=f'{LatestAssembly}_genomic.fna.gz'
    #download genome
    call(f'wget --output-document=../../data/genomes/gDNA.fna.gz {URL+filename}',shell=True)
    #decompress genome
    call(f'gzip -d ../../data/genomes/gDNA.fna.gz',shell=True)
    #get mitochondrial id
    mtID=run(
        'egrep \"Homo sapiens mitochondrion\" ../../data/genomes/mitochondrion.1.1.genomic.fna',
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
    #generate LASTAL db
    call('lastdb ../../data/genomes/db ../../data/genomes/gDNA.fna', shell=True)
    #align gDNA with dmtDNA
    call('lastal -r1 -q1 -a7 -b1 ../../data/genomes/db dmtDNA.fna > dmtDNA.afa', shell=True)
    #filter significant alignments



#apply function
proba=organisms.apply(process_DNA)