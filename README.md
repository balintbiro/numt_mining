# NUMT mining pipeline
Pipeline for mining mitochondrial sequences in the nuclear genome.

# Upon usage, please cite
 ``` bib
 @article{biro2024mitochondrial,
  title={Mitochondrial genome plasticity of mammalian species},
  author={Bir{\'o}, B{\'a}lint and G{\'a}l, Zolt{\'a}n and Fekete, Zs{\'o}fia and Klecska, Eszter and Hoffmann, Orsolya Ivett},
  journal={BMC Genomics},
  volume={25},
  number={1},
  pages={1--14},
  year={2024},
  publisher={BioMed Central}
}
 ```

# Abstract
There is an ongoing process in which mitochondrial sequences are being integrated into the nuclear genome. The importance of these sequences has already been revealed in cancer biology, forensic, phylogenetic studies and in the evolution of the eukaryotic genetic information. Human and numerous model organisms’ genomes were described from those sequences point of view. Furthermore, recent studies were published on the patterns of these nuclear localised mitochondrial sequences in different taxa.

However, the results of the previously released studies are difficult to compare due to the lack of standardised methods and/or using few numbers of genomes. Therefore, in this paper our primary goal is to establish a uniform mining pipeline to explore these nuclear localised mitochondrial sequences.

Our results show that the frequency of several repetitive elements is higher in the flanking regions of these sequences than expected. A machine learning model reveals that the flanking regions' repetitive elements and different structural characteristics are highly influential during the integration process.

In this paper, we introduce a general mining pipeline for all mammalian genomes. The workflow is publicly available and is believed to serve as a validated baseline for future research in this field. We confirm the widespread opinion, on - as to our current knowledge - the largest dataset, that structural circumstances and events corresponding to repetitive elements are highly significant. An accurate model has also been trained to predict these sequences and their corresponding flanking regions

![graphical_abstract](/results/fig6.png)

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
