#loading the required packages
library(seqinr)
library(phangorn)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

#read in the previously prepared sequences
sequences <- readDNAStringSet('Documents/Projects/numt_mining/cob_genes_and_numts_small_header.txt',
                              format='fasta')
#reorient the sequences
sequences <- OrientNucleotides(sequences)

#align the sequences
aligned_sequences <- AlignSeqs(sequences)

#write the aligned sequences into a file
writeXStringSet(aligned_sequences,
                file='Documents/Projects/numt_mining/aligned_cobs_small_header.fasta')

#read the alignment
dna <- read.alignment('Documents/Projects/numt_mining/aligned_cobs_small_header.fasta',
                      format = 'fasta')
#creating distance matrix
distance_matrix <- dist.alignment(dna)

treeupgma <- upgma(distance_matrix)
treenj <- nj(dm)

treepars <- optim.parsimony(treeupgma,as.phyDat(dna))
treeratchet <- pratchet(as.phyDat(dna),trace=0)
parsimony(c(treepars, treeratchet), as.phyDat(dna))
treeratchet <- acctran(treeratchet, as.phyDat(dna))

ggtree(root(treenj,outgroup='rat_CYTB'))+
  geom_tiplab(hjust = -.5, size=3.5, align = TRUE)+
  xlim(0,1)

ggtree(root(treeupgma,outgroup='rat_CYTB'))+
  geom_tiplab(hjust = -.5, size=3.5, align = TRUE)+
  xlim(0,1)

ggtree(root(treeratchet,outgroup='rat_CYTB'))+
  geom_tiplab(hjust = -.5, size=3.5, align = TRUE)+
  xlim(0,180)

ggtree(root(treepars,outgroup='rat_CYTB'))+
  geom_tiplab(hjust = -.5, size=3.5, align = TRUE)+
  xlim(0,10)
