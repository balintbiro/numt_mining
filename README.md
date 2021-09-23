# numt_mining

Pipeline for mining mitochondrial sequences in the nuclear genome.
Multiple sequence alignments are performed with LASTAL.
The LINUX commands are embedded into Python scripts with the subprocess module.
Nuclear and mitochondrial genomes are acquired from Ensembl genome browser (http://ftp.ensembl.org/pub/release-104/fasta/).
The initial folder structure consists of three folders:

- code
- data
- results

By default, the data folder contains a plain text file called organism_names.txt.
In this file the latin names of the organisms should be provided with capital first letter, '_' sign in between and separated by ',' commas
(e.g. Homo_sapiens,Mus_musculus). 

Used softwares are the followings (version number is provided if available):

- LASTAL (1219)
- Python (3.7.10)

Used Python modules are the followings (version number is provided if available):

- os
- gzip
- shutil
- pandas (1.2.4)
- subprocess
- urllib.request (3.7)

Programs:

data_acquisition.py

- Download nuclear (gDNA) and mitochondrial (mtDNA) sequences from Ensembl.
- If no gDNA and/or mtDNA available, a subfolder called problem_reports with the organism name will be created inside the results folder.
Separated txt files will be written in the problem reports folder with the corresponding organism names.

alignment.py

- Create mtDNA permutations (reversed and doubled).
- Build the required LASTAL database.
- Do the alignment between gDNA and reversed mtDNA.
- Get the e-value threshold.
- Align double mtDNA and gDNA.
- Filter the alignment based on the precomputed e-value threshold.