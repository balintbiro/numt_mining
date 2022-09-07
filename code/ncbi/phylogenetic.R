#loading dependencies
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
library(phangorn)

#set working directory
getwd()
setwd('Documents/Projects/numt_mining/code/ncbi/')

###############################################################Neighbor Joining Tree

#read in alignment (previously generated with clustal omega)
alignment <- read.alignment(
  '../../data/aligned_mtDNAs.fa',
  format="fasta"
)

#calculate distance matrix
distance_matrix <- dist.alignment(
  alignment,
  matrix="similarity"
)

#create the actual tree
nj_tree <- nj(distance_matrix)

#visualize nj tree
ggtree(
  nj_tree,
  layout='fan',
  branch.length='none',
  size=0.1,
  color='black',
  alpha=1,
  ,open.angle=90
)+
  geom_tippoint(
    size=.15,
    color='black',
    alpha=.5
  )+theme(
    plot.background = element_rect(
      fill='transparent',
      color=NA
    ),
    panel.background = element_rect(
      fill='transparent',
      color=NA
    )
  )+geom_tiplab(
    hjust=.4
  )#+geom_highlight(node=1,fill='purple',alpha=.9)

#save nj tree
ggsave(
  '../../results/nj_tree.png',
  plot=last_plot(),
  dpi=400,units='cm',
  width=14,
  height=14,
  bg='transparent'
)

###############################################################UPGMA Tree

#read alignment
alignment <- read.phyDat(
  '../../data/aligned_mtDNAs.fa',
  format='fasta'
)

#create distance matrix
distance_matrix <- dist.ml(alignment)

#create the tree itself
upgma_tree <- upgma(distance_matrix)

#visualize tree
ggtree(
  upgma_tree,
  layout='fan',
  branch.length='none',
  size=1.8,
  color='black',
  alpha=1
)+#,open.angle=180 can be used together with fan layout
  geom_tippoint(
    size=.5,
    color='black',
    alpha=.5
  )+geom_highlight(
    node=1,
    fill='purple',
    alpha=.9
  )
