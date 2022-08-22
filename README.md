# numt mining pipeline
Pipeline for mining mitochondrial sequences in the nuclear genome.

Setting up the environment with conda package manager
---
```bash
mkdir codes data results
conda install -c bioconda last=1219
conda install -c bioconda samtools=1.6
conda install -c bioconda repeatmasker=4.1.2
```
In the codes relative paths are used so in case of forking the repository it is necessary to have the same folder structure as described above.

Nuclear genome sequence files are from
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/ , while the files (mitochondrion.1.1.genomic.fna.gz, mitochondrion.2.1.genomic.fna.gz) containing the mitochondrial sequences are from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/ . Tax  information is from https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt and phylogenetic ranks are from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip / rankedlineage.dmp.