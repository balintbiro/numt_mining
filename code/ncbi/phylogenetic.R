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
bs_tree <- bootstrap.phyDat(as.phyDat(alignment),bs=500,fun)

#visualize the consensus tree
tree <- ggtree(
  consensus(bs_tree),#the consensus tree
  layout='circular',#style,other options are daylight, equal_angle, fan, circular, roundrect, ellipse, slanted, unrooted
  branch.length='none',#no branch length
  size=0.4,#width of the branches
  color='black',#color
  alpha=1,#transparency
  open.angle=90#where to have a gap
)+geom_highlight(#higlight the node with box
  node=156,#node number is for Carnivores
  fill='#1f77b4ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75,
  type="rect"
)+geom_highlight(#higlight the node with box
  node=182,#node number is for Primates
  fill='#aec7e8ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=247,#node number is for Primates
  fill='#aec7e8ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=198,#node number is for Rodentia
  fill='#2ca02cff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=248,#node number is for Rodentia
  fill='#2ca02cff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=244,#node number is for Rodentia
  fill='#2ca02cff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=252,#node number is for Rodentia
  fill='#2ca02cff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=228,#node number is for Chiroptera
  fill='#ff7f0eff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=235,#node number is for Chiroptera
  fill='#ff7f0eff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=250,#node number is for Chiroptera
  fill='#ff7f0eff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=208,#node number is for Artiodactyla
  fill='#ffbb78ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=216,#node number is for Artiodactyla
  fill='#ffbb78ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=232,#node number is for Artiodactyla
  fill='#ffbb78ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=226,#node number is for Didelphimorpha
  fill='#e377c2ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=243,#node number is for Perissodactyla
  fill='#8c564bff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=240,#node number is for Lagomopha
  fill='#c7c7c7ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=227,#node number is for Monotremata
  fill='#bcbd22ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=246,#node number is for Proboscidea
  fill='#f7b6d2ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=238,#node number is for Eulipotyphla
  fill='#d62728ff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)+geom_highlight(#higlight the node with box
  node=251,#node number is for Pholidota
  fill='#7f7f7fff',#color to fill with;rgb and hex are also accepted; alpha is also possible to be specified
  alpha=.75
)




#modify background
theme(
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
  size=4#size of the text
)+geom_text2(
  aes(subset=!isTip,label=node)
)

#add nodelabel as a stick
tree+geom_cladelabel(
  node=156,#node number
  label='test label',#text for highlighting
  color='green',#color to fill with
  fontsize=8,#size of the text
  angle=45,
  barsize=7#thickness of bar
)

ggtree(
  consensus(bs_tree),#the consensus tree
  layout='circular',#style,other options are daylight, equal_angle, fan, circular, roundrect, ellipse, slanted, unrooted
  branch.length='none',#no branch length
  size=0.4,#width of the branches
  color='black',#color
  alpha=1,#transparency
  open.angle=90#where to have a gap
)+geom_cladelabel(#higlight the node with box
  node=156,#node number is for Carnivores
  color='#1f77b4ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=182,#node number is for Primates
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=247,#node number is for Primates
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=198,#node number is for Rodentia
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=248,#node number is for Rodentia
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=244,#node number is for Rodentia
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=252,#node number is for Rodentia
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=228,#node number is for Chiroptera
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=235,#node number is for Chiroptera
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=250,#node number is for Chiroptera
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=208,#node number is for Artiodactyla
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=216,#node number is for Artiodactyla
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=232,#node number is for Artiodactyla
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=226,#node number is for Didelphimorpha
  color='#e377c2ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=243,#node number is for Perissodactyla
  color='#8c564bff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=240,#node number is for Lagomopha
  color='#c7c7c7ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=227,#node number is for Monotremata
  color='#bcbd22ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=246,#node number is for Proboscidea
  color='#f7b6d2ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=238,#node number is for Eulipotyphla
  color='#d62728ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=251,#node number is for Pholidota
  color='#7f7f7fff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)

#save nj tree
ggsave(
  '../../results/bs_nj_tree.png',
  plot=last_plot(),
  dpi=400,
  units='cm',
  width=14,
  height=14,
  bg='transparent'
)
