# numt_mining

Pipeline for mining mitochondrial sequences in the nuclear genome. Multiple sequence alignments are done with LASTAL. The LINUX commands for LASTAL are embedded into Python scripts. Nuclear and mitochondrial genomes are acquired from Ensembl genome browser (http://ftp.ensembl.org/pub/release-104/fasta/).

Used Python modules are the followings (version number is provided if available):

- os
- gzip
- shutil
- pandas (1.2.4)
- subprocess
- urllib.request (3.7)
