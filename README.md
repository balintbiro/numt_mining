# numt_mining

Pipeline for mining mitochondrial sequences in the nuclear genome.
Multiple sequence alignments are performed with LASTAL.
The LINUX commands are embedded into Python scripts with the subprocess module.
Nuclear and mitochondrial genomes are acquired from Ensembl genome browser (http://ftp.ensembl.org/pub/release-104/fasta/). Mitochondrial fasta sequences are also available at https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/ link!
The initial folder structure consists of three folders:

- code
- data
- results

By default, the data folder contains a plain text file called organism_names.txt.
In this file the latin names of the organisms should be provided with capital first letter, '_' sign in between and separated by ',' commas
(e.g. Homo_sapiens,Mus_musculus). 

Used softwares are the followings (version number is provided if available):

- GNU Wget (1.14)
- grep
- LASTAL (1219)
	- https://github.com/UCSantaCruzComputationalGenomicsLab/last/blob/master/doc/lastal.txt
	- http://last.cbrc.jp/doc/lastal.txt
- Python (3.7.10)

Used Python modules are the followings (version number is provided if available):

- os
- gzip
- shutil
- pandas (1.2.4)
- filecmp
- subprocess
- urllib.request (3.7)

Programs:

data_acquisition.py

- Download nuclear (gDNA) and mitochondrial (mtDNA) sequences from Ensembl.
- If no gDNA and/or mtDNA available, a subfolder called problem_reports with the organism name will be created inside the results folder.
Separated txt files will be written in the problem reports folder with the corresponding organism names.

duplicate_mt.py

- Get the mitochondrial sequence from mt.fa and create a duplicated version (d_mt) of it
- Write the d_mt.fa file

alignment.py

- Create mtDNA permutations (reversed and doubled).
- Build the required LASTAL database.
- Align double mtDNA and gDNA.


Nuclear genome informations (size and quality) can be downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/ as a .txt file called eukaryotes.txt, while mitochondrial genome informations can be obtained from https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/ as a .csv file called organelles.csv. Phylogenetic information is gained from https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip included nodes.dmp file. Mammalian genomes are also available at https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/.