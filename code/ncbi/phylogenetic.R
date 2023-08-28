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

#######new
distace_matrix <- dist.ml(as.phyDat(alignment))
treeNJ <- NJ(distance_matrix)
fit <- pml(treeNJ,data=as.phyDat(alignment))
fitJC <- optim.pml(fit,rearrangement = 'NNI')
fitGTR <- optim.pml(fitJC,model='GTR',optInv = TRUE,optGamma = TRUE,
                    rearrangement = 'NNI',control=pml.control(trace=0))
bs <- bootstrap.pml(fitGTR, bs=100, optNni=TRUE,
                    control = pml.control(trace = 0))

tree<-ggtree(consensus(bs),layout='circular')

write.tree(consensus(bs),file='../../results/gtr_tree.nex')

trial_tree<-read.tree('../../results/gtr_tree_trial.nex')

ggtree(trial_tree,layout='circular',size=0.4,open.angle = 280)+#geom_tiplab()+geom_text2(aes(subset=!isTip,label=node))+
  geom_cladelabel(node=236,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=253,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=254,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=256,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=257,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=163,color='#aec7e8ff',barsize=7,label='')+#Primates
  geom_cladelabel(node=197,color='#aec7e8ff',barsize=7,label='')+#Primates
  geom_cladelabel(node=178,color='#2ca02cff',barsize=7,label='')+#Rodentia
  geom_cladelabel(node=198,color='#2ca02cff',barsize=7,label='')+#Rodentia
  geom_cladelabel(node=194,color='#2ca02cff',barsize=7,label='')+#Rodentia
  geom_cladelabel(node=204,color='#ffbb78ff',barsize=7,label='')+#Artiodactyla
  geom_cladelabel(node=225,color='#ff7f0eff',barsize=7,label='')+#Chiroptera
  geom_cladelabel(node=235,color='#ff7f0eff',barsize=7,label='')+#Chiroptera
  geom_cladelabel(node=232,color='#ff7f0eff',barsize=7,label='')+#Chiroptera
  geom_cladelabel(node=199,color='#c7c7c7ff',barsize=7,label='')+#Lagomorpha
  geom_cladelabel(node=202,color='#f7b6d2ff',barsize=7,label='')+#Proboscidea
  geom_cladelabel(node=190,color='#e377c2ff',barsize=7,label='')+#Didelphimorpha
  geom_cladelabel(node=192,color='#e377c2ff',barsize=7,label='')+#Didelphimorpha
  geom_cladelabel(node=230,color='#8c564bff',barsize=7,label='')+#Perissodactyla
  geom_cladelabel(node=228,color='#d62728ff',barsize=7,label='')+#Eulipotyphla
  geom_cladelabel(node=234,color='#7f7f7fff',barsize=7,label='')+#Pholidota
  geom_cladelabel(node=193,color='#bcbd22ff',barsize=7,label='')+#Monotremata
  geom_cladelabel(node=178,color='#2ca02cff',barsize=7,label='')#Rodentia




#save nj tree
ggsave(
  '../../results/gtr_tree.png',
  plot=last_plot(),
  dpi=400,
  units='cm',
  width=14,
  height=14,
  bg='transparent'
)
