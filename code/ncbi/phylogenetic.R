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

ggtree(trial_tree,layout='circular',size=0.4,open.angle = 280)+geom_tiplab()+geom_text2(aes(subset=!isTip,label=node))+
  geom_cladelabel(node=236,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=253,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=254,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=256,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=257,color='#1f77b4ff',barsize=7,label='')+#Carnivores
  geom_cladelabel(node=163,color='#aec7e8ff',barsize=7,label='')+#Primates
  geom_cladelabel(node=197,color='#aec7e8ff',barsize=7,label='')+#Primates
  geom_cladelabel(node=198,color='#2ca02cff',barsize=7,label='')+#Rodentia
  geom_cladelabel(node=194,color='#2ca02cff',barsize=7,label='')+#Rodentia
  geom_cladelabel(node=204,color='#ffbb78ff',barsize=7,label='')+#Artiodactyla
  geom_cladelabel(node=225,color='#ff7f0eff',barsize=7,label='')+#Chiroptera
  geom_cladelabel(node=232,color='#ff7f0eff',barsize=7,label='')+#Chiroptera
  geom_cladelabel(node=178,color='#2ca02cff',barsize=7,label='')#Rodentia



++geom_cladelabel(#higlight the node with box
  node=217,#node number is for Rodentia=['castor_canadensis', 'cavia_porcellus', 'chinchilla_lanigera', 'cricetulus_griseus', 'fukomys_damarensis', 'grammomys_surdaster', 'heterocephalus_glaber', 'ictidomys_tridecemlineatus', 'jaculus_jaculus', 'marmota_flaviventris', 'mastomys_coucha', 'meriones_unguiculatus', 'mesocricetus_auratus', 'microtus_ochrogaster', 'mus_caroli', 'mus_musculus', 'mus_pahari', 'myodes_glareolus', 'nannospalax_galili', 'octodon_degus', 'peromyscus_leucopus', 'rattus_norvegicus', 'rattus_rattus', 'urocitellus_parryii']
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=235,#node number is for Rodentia=['castor_canadensis', 'cavia_porcellus', 'chinchilla_lanigera', 'cricetulus_griseus', 'fukomys_damarensis', 'grammomys_surdaster', 'heterocephalus_glaber', 'ictidomys_tridecemlineatus', 'jaculus_jaculus', 'marmota_flaviventris', 'mastomys_coucha', 'meriones_unguiculatus', 'mesocricetus_auratus', 'microtus_ochrogaster', 'mus_caroli', 'mus_musculus', 'mus_pahari', 'myodes_glareolus', 'nannospalax_galili', 'octodon_degus', 'peromyscus_leucopus', 'rattus_norvegicus', 'rattus_rattus', 'urocitellus_parryii']
  color='#2ca02cff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=194,#node number is for Artiodactyla=['balaenoptera_musculus', 'bos_indicus', 'bos_mutus', 'bos_taurus', 'bubalus_bubalis', 'camelus_bactrianus', 'camelus_dromedarius', 'camelus_ferus', 'capra_hircus', 'cervus_canadensis', 'cervus_elaphus', 'delphinapterus_leucas', 'globicephala_melas', 'lagenorhynchus_obliquidens', 'lipotes_vexillifer', 'monodon_monoceros', 'odocoileus_virginianus', 'orcinus_orca', 'ovis_aries', 'phacochoerus_africanus', 'phocoena_sinus', 'physeter_catodon', 'sus_scrofa', 'tursiops_truncatus', 'vicugna_pacos']
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=238,#node number is for Artiodactyla=['balaenoptera_musculus', 'bos_indicus', 'bos_mutus', 'bos_taurus', 'bubalus_bubalis', 'camelus_bactrianus', 'camelus_dromedarius', 'camelus_ferus', 'capra_hircus', 'cervus_canadensis', 'cervus_elaphus', 'delphinapterus_leucas', 'globicephala_melas', 'lagenorhynchus_obliquidens', 'lipotes_vexillifer', 'monodon_monoceros', 'odocoileus_virginianus', 'orcinus_orca', 'ovis_aries', 'phacochoerus_africanus', 'phocoena_sinus', 'physeter_catodon', 'sus_scrofa', 'tursiops_truncatus', 'vicugna_pacos']
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=201,#node number is for Artiodactyla=['balaenoptera_musculus', 'bos_indicus', 'bos_mutus', 'bos_taurus', 'bubalus_bubalis', 'camelus_bactrianus', 'camelus_dromedarius', 'camelus_ferus', 'capra_hircus', 'cervus_canadensis', 'cervus_elaphus', 'delphinapterus_leucas', 'globicephala_melas', 'lagenorhynchus_obliquidens', 'lipotes_vexillifer', 'monodon_monoceros', 'odocoileus_virginianus', 'orcinus_orca', 'ovis_aries', 'phacochoerus_africanus', 'phocoena_sinus', 'physeter_catodon', 'sus_scrofa', 'tursiops_truncatus', 'vicugna_pacos']
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=213,#node number is for Chiroptera=['artibeus_jamaicensis', 'desmodus_rotundus', 'hipposideros_armiger', 'myotis_brandtii', 'myotis_davidii', 'myotis_lucifugus', 'myotis_myotis', 'pteropus_alecto', 'pteropus_vampyrus', 'rhinolophus_ferrumequinum', 'rousettus_aegyptiacus']
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=223,#node number is for Chiroptera=['artibeus_jamaicensis', 'desmodus_rotundus', 'hipposideros_armiger', 'myotis_brandtii', 'myotis_davidii', 'myotis_lucifugus', 'myotis_myotis', 'pteropus_alecto', 'pteropus_vampyrus', 'rhinolophus_ferrumequinum', 'rousettus_aegyptiacus']
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=237,#node number is for Chiroptera=['artibeus_jamaicensis', 'desmodus_rotundus', 'hipposideros_armiger', 'myotis_brandtii', 'myotis_davidii', 'myotis_lucifugus', 'myotis_myotis', 'pteropus_alecto', 'pteropus_vampyrus', 'rhinolophus_ferrumequinum', 'rousettus_aegyptiacus']
  color='#ff7f0eff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=220,#node number is for Artiodactyla=['balaenoptera_musculus', 'bos_indicus', 'bos_mutus', 'bos_taurus', 'bubalus_bubalis', 'camelus_bactrianus', 'camelus_dromedarius', 'camelus_ferus', 'capra_hircus', 'cervus_canadensis', 'cervus_elaphus', 'delphinapterus_leucas', 'globicephala_melas', 'lagenorhynchus_obliquidens', 'lipotes_vexillifer', 'monodon_monoceros', 'odocoileus_virginianus', 'orcinus_orca', 'ovis_aries', 'phacochoerus_africanus', 'phocoena_sinus', 'physeter_catodon', 'sus_scrofa', 'tursiops_truncatus', 'vicugna_pacos']
  color='#ffbb78ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=234,#node number is for Primates=['aotus_nancymaae', 'callithrix_jacchus', 'carlito_syrichta', 'cercocebus_atys', 'chlorocebus_sabaeus', 'gorilla_gorilla', 'lemur_catta', 'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'mandrillus_leucophaeus', 'microcebus_murinus', 'nomascus_leucogenys', 'pan_paniscus', 'pan_troglodytes', 'papio_anubis', 'pongo_abelii', 'propithecus_coquereli', 'rhinopithecus_bieti', 'rhinopithecus_roxellana', 'sapajus_apella', 'theropithecus_gelada', 'trachypithecus_francoisi']
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=233,#node number is for Proboscidea=['elephas_maximus', 'loxodonta_africana']
  color='#f7b6d2ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=231,#node number is for Perissodactyla=['equus_asinus', 'equus_caballus']
  color='#8c564bff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=228,#node number is for Lagomopha=['ochotona_curzoniae', 'ochotona_princeps', 'oryctolagus_cuniculus']
  color='#c7c7c7ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=226,#node number is for Eulipotyphla=['condylura_cristata', 'erinaceus_europaeus', 'sorex_araneus', 'talpa_occidentalis']
  color='#d62728ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=239,#node number is for Pholidota=['manis_javanica', 'manis_pentadactyla']
  color='#7f7f7fff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=211,#node number is for Didelphimorpha=['gracilinanus_agilis', 'monodelphis_domestica']
  color='#e377c2ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=212,#node number is for Monotremata=['ornithorhynchus_anatinus', 'tachyglossus_aculeatus']
  color='#bcbd22ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=240,#node number is for Primates=['aotus_nancymaae', 'callithrix_jacchus', 'carlito_syrichta', 'cercocebus_atys', 'chlorocebus_sabaeus', 'gorilla_gorilla', 'lemur_catta', 'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'mandrillus_leucophaeus', 'microcebus_murinus', 'nomascus_leucogenys', 'pan_paniscus', 'pan_troglodytes', 'papio_anubis', 'pongo_abelii', 'propithecus_coquereli', 'rhinopithecus_bieti', 'rhinopithecus_roxellana', 'sapajus_apella', 'theropithecus_gelada', 'trachypithecus_francoisi']
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=241,#node number is for Primates=['aotus_nancymaae', 'callithrix_jacchus', 'carlito_syrichta', 'cercocebus_atys', 'chlorocebus_sabaeus', 'gorilla_gorilla', 'lemur_catta', 'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'mandrillus_leucophaeus', 'microcebus_murinus', 'nomascus_leucogenys', 'pan_paniscus', 'pan_troglodytes', 'papio_anubis', 'pongo_abelii', 'propithecus_coquereli', 'rhinopithecus_bieti', 'rhinopithecus_roxellana', 'sapajus_apella', 'theropithecus_gelada', 'trachypithecus_francoisi']
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)+geom_cladelabel(#higlight the node with box
  node=250,#node number is for Primates=['aotus_nancymaae', 'callithrix_jacchus', 'carlito_syrichta', 'cercocebus_atys', 'chlorocebus_sabaeus', 'gorilla_gorilla', 'lemur_catta', 'macaca_fascicularis', 'macaca_mulatta', 'macaca_nemestrina', 'mandrillus_leucophaeus', 'microcebus_murinus', 'nomascus_leucogenys', 'pan_paniscus', 'pan_troglodytes', 'papio_anubis', 'pongo_abelii', 'propithecus_coquereli', 'rhinopithecus_bieti', 'rhinopithecus_roxellana', 'sapajus_apella', 'theropithecus_gelada', 'trachypithecus_francoisi']
  color='#aec7e8ff',#color to color with;rgb and hex are also accepted; alpha is also possible to be specified
  barsize=7,
  label=''
)


#save nj tree
ggsave(
  '../../results/bs_nj_tree.eps',
  plot=last_plot(),
  dpi=400,
  units='cm',
  width=14,
  height=14,
  bg='transparent'
)
