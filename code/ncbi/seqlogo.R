#installing ggseqlogo package
#install.packages('ggseqlogo')

#load dependencies
library(ggseqlogo)
library(ggplot2)
library(seqinr)

setwd('Documents/Projects/numt_mining/code/ncbi/')

system(
  command='ls ../'
)

sequences <- read.fasta(
  file='../../data/aligned_seqlogos.fa',
  seqonly=TRUE
)
ggseqlogo(
  unlist(sequences),
  seq_type='dna'
)
