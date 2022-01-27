#loading the required packages
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)

#read in the previously prepared sequences
sequences <- readDNAStringSet('Documents/Projects/numt_mining/cob_genes_and_numts.txt',
                              format='fasta')
#reorient the sequences
sequences <- OrientNucleotides(sequences)

#align the sequences
aligned_sequences <- AlignSeqs(sequences)

#open up the aligned sequences in a browser tab
#BrowseSeqs(aligned_sequences,
#           highlight = 0)

#write the aligned sequences into a file
writeXStringSet(aligned_sequences,
                file='Documents/Projects/numt_mining/aligned_cobs.fasta')

#read the alignment
dna <- read.alignment('Documents/Projects/numt_mining/aligned_cobs.fasta',
                      format = 'fasta')
#creating distance matrix
distance_matrix <- dist.alignment(dna,
                    matrix = 'similarity')

#using njs() function because of the missing values, it turned our distance matrix ito a phylo object
tree <- nj(distance_matrix)
class(tree)

options(ignore.negative.edge=TRUE)

fig <- ggtree(root(tree,outgroup = 'rat_cob'))+
  geom_tiplab(hjust = -.5, size=3.5, align = TRUE) +
  xlim(0,1.35)+
  geom_hilight(node=1, fill="purple", alpha = 0.2) +
  geom_cladelabel(node=35, label="outgroup", 
                  color="blue", offset=.57, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5, hjust=-.1)+
  geom_cladelabel(node=52, label="mt-CYTB", 
                  color="dark green", offset=.57, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5, hjust=-.1)+
  geom_cladelabel(node=53, label="CYTB-numt", 
                  color="brown", offset=.57, barsize = 2,
                  fontsize = 5, align=TRUE, alpha = 0.5, hjust=-.1)
ggsave('Documents/Projects/numt_mining/results//pyholgenetic_cobs_numts.png',
       plot=last_plot(),
       , width = 7, height = 4)

