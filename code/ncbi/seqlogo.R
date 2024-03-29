#installing ggseqlogo package
#install.packages('ggseqlogo')

#load dependencies
library(ggseqlogo)
library(ggplot2)
library(seqinr)

#set working directory
setwd('Documents/Projects/numt_mining/code/ncbi/')

#read inputs
seqlogo_input <- read.csv(
  file='../../data/seqlogo_inputs.csv'
)

#check if there is a directory for the seqlogos
if (file.exists('../../results/seqlogos/')==FALSE){
  dir.create('../../results/seqlogos/')
}

#function for visualizing
order_seqlogo_vis <- function(order){
  fil <- seqlogo_input[,'order']==order
  seqlogo <- ggseqlogo(
    sequences[fil,]$seq,
    seq_type='dna',
    method='bits'
  )+xlab('Positions')+
    theme(axis.title.x=element_text(size=30),
          axis.title.y=element_text(size=30))
  ggsave(
    sprintf('../../results/seqlogos/%s.png',order),
    device=NULL,
    plot=seqlogo,
    bg='white',
    dpi=300
  )
}

order_seqlogo_vis('Monotremata')

repeat_seqlogo_vis <- function(repeatclass){
  fil <- seqlogo_input[,'u_1st_repeatclass']==repeatclass
  seqlogo <- ggseqlogo(
    sequences[fil,]$seq,
    seq_type='dna',
    method='bits'
  )+xlab('Positions')+
    theme(axis.title.x=element_text(size=30),
          axis.title.y=element_text(size=30))
  ggsave(
    sprintf('../../results/seqlogos/%s.png',gsub('/','_',repeatclass)),
    device=NULL,
    plot=seqlogo,
    bg='white',
    dpi=300
  )
}

repeat_seqlogo_vis('SINE/B4')
