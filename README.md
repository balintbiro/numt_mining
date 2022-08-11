# numt mining pipeline
Pipeline for mining mitochondrial sequences in the nuclear genome.
The initial folder structure consists of three folders:

- code
- data
- results

Setting up the environment
---
```bash
mkdir codes data results
conda install -c bioconda last
conda install -c bioconda repeatmasker
```

Multiple sequence alignments are performed with LASTAL.
The LINUX commands are embedded into Python scripts with the subprocess module.
Nuclear genome sequence files were acquired from
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/ , while the files (mitochondrion.1.1.genomic.fna.gz, mitochondrion.2.1.genomic.fna.gz) containing the mitochondrial sequences were acquired from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/ .