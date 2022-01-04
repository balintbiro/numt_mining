#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#load the package
library(biomaRt)

#check the avaialable databases
#listEnsembl()
#check the avaiable datasets
#listDatasets(ensembl)
ensembl <- useEnsembl(biomart='mouse_strains',#select the desired database
                      dataset='mm129s1svimj_gene_ensembl',#select the desired dataset
                      mirror='asia')#select the mirror, options -->useast,uswest,asia,www

#check the available filters and attributes
#listFilters(ensembl)
#listAttributes(ensembl)
genes<-getBM(attributes=c('external_gene_name','description','uniprotswissprot'),
      filters=c('chromosome_name','start','end'),
      values=list(chromosome='4',start='71110168',end='71110178'),
      mart=ensembl)
