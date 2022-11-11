#get node numbers https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/ggtree/inst/doc/treeManipulation.html#internal-node-number
#clade annotations https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html
#possible bootstrapping for densitree fun <- function(x) nj(dist.ml(x,model='JC69')); bootstrap.phyDat(alignment,  fun, bs=100)
#visualize bootstrap tree-plotBS()
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
tree <- ggtree(
  nj_tree,
  layout='fan',#other options are circular, roundrect, ellipse, slanted, unrooted
  branch.length='none',
  size=0.1,
  color='black',
  alpha=1,
  ,open.angle=90
)

#add points to the ends of branches
tree+geom_tippoint(
    size=.15,
    color='black',
    alpha=.5
  )

#modify background
tree+theme(
    plot.background = element_rect(
      fill='transparent',
      color=NA
    ),
    panel.background = element_rect(
      fill='transparent',
      color=NA
    )
  )

#add organism names-it is going to be messy since we have too much names
tree+geom_tiplab(
    hjust=.4
  )

#save nj tree
ggsave(
  '../../results/nj_tree.png',
  plot=last_plot(),
  dpi=400,units='cm',
  width=14,
  height=14,
  bg='transparent'
)
###############################################################add node numbers
#visualize nj tree
tree <- ggtree(
  nj_tree,
  layout='fan',#other options are circular, roundrect, ellipse, slanted, unrooted
  branch.length='none',
  size=0.1,
  color='black',
  alpha=1,
  ,open.angle=90
)

#add points to the ends of branches
tree+geom_tippoint(
  size=.15,
  color='black',
  alpha=.5
)

#modify background
tree+theme(
  plot.background = element_rect(
    fill='transparent',
    color=NA
  ),
  panel.background = element_rect(
    fill='transparent',
    color=NA
  )
)

#add organism names-it is going to be messy since we have too much names
tree+geom_tiplab(
  hjust=.4
)

#add node numbers
tree+geom_text2(
  aes(subset=!isTip,label=node)
)

#add nodelabel as a stick
tree+geom_cladelabel(
  node=158,label='test label',color='green',fontsize=8,angle=45,barsize=3
)

#highlight node with box
tree+geom_highlight(
  node=161, fill='steelblue',alpha=.3
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
