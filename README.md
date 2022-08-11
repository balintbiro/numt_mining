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

Nuclear genome sequence files are acquired from
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/ , while the files (mitochondrion.1.1.genomic.fna.gz, mitochondrion.2.1.genomic.fna.gz) containing the mitochondrial sequences are acquired from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/ .

Used LINUX commands
---
```bash
wget
gzip
```

Used Python packages
---
```python
subprocess
os
```