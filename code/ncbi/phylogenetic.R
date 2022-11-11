#get node numbers https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/ggtree/inst/doc/treeManipulation.html#internal-node-number
#clade annotations https://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html
#possible bootstrapping for densitree fun <- function(x) nj(dist.ml(x,model='JC69')); bootstrap.phyDat(alignment,  fun, bs=100)
#visualize bootstrap tree-plotBS
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

#function for bootstrapping a neighbor joining tree with JC method
fun <- function(x) nj(dist.ml(x,model='JC69'))

#actual bootstrapping 100 times
bs_tree <- bootstrap.phyDat(as.phyDat(alignment),bs=100,fun)

#visualize the consensus tree
tree <- ggtree(
  consensus(bs_tree),#the consensus tree
  layout='circular',#style,other options are daylight, fan, circular, roundrect, ellipse, slanted, unrooted
  branch.length='none',#no branch length
  size=0.1,#width of the branches
  color='black',#color
  alpha=1,#transparency
  open.angle=90#where to have a gap
)

#add points to the ends of branches
tree+geom_tippoint(
  size=.4,
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
  hjust=-.1,#distance from the tip
  size=2#size of the text
)

#add node numbers
tree+geom_text2(
  aes(subset=!isTip,label=node)
)

#add nodelabel as a stick
tree+geom_cladelabel(
  node=214,#node number
  label='test label',#text for highlighting
  color='green',#color to fill with
  fontsize=8,#size of the text
  angle=45,
  barsize=3#thickness of bar
)

#highlight node with box
tree+geom_highlight(
  node=179,#node number
  fill='#000080',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
)

#save nj tree
ggsave(
  '../../results/nj_tree.png',
  plot=last_plot(),
  dpi=400,
  units='cm',
  width=14,
  height=14,
  bg='transparent'
)