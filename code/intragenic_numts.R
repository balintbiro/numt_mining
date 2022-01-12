
#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#load the package
library(biomaRt)

ensembl_db<-useEnsembl(biomart='mouse_strains',#select the desired database, all available listEnsembl()
                    dataset='mmaj_gene_ensembl',#select the desired dataset, all available listDatasets(ensembl)
                    mirror='useast')#select the mirror, options -->useast,uswest,asia,www

genomic_regions<-scan('Z:/balint/numt/mice_strains_numt/genomic_regions/mmaj_gene_ensembl.txt',character())

ensembl_query<-function(genomic_region){
  ensembl<-getBM(attributes=c('external_gene_name'),
                 filters=c('chromosomal_region'),#all available listFilters(ensembl)
                 values=list(genomic_region),#all available listAttributes(ensembl)
                 mart=ensembl_db)
  return(ensembl)
}

genes<-sapply(genomic_regions,ensembl_query)
genes<-unlist(genes)
df<-stack(genes)