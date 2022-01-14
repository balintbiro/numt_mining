
#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#install.packages('BiocManager',quietly=TRUE)
BiocManager::install('biomaRt')
#load the package
library(biomaRt)

current_strain <- 'mmwsbeij'

ensembl_db<-useEnsembl(biomart='mouse_strains',#select the desired database, all available listEnsembl()
                    dataset=sprintf('%s_gene_ensembl',current_strain),#select the desired dataset, all available listDatasets(ensembl)
                    mirror='uswest')#select the mirror, options -->useast,uswest,asia,www

#get the strain specific genomic regions
genomic_regions<-scan(sprintf('/Volumes/motilin/balint/numt/mice_strains_numt/genomic_regions/%s_gene_ensembl.txt',current_strain),character())

#function for querying Ensembl
ensembl_query<-function(genomic_region){
  ensembl<-getBM(attributes=c('external_gene_name'),
                 filters=c('chromosomal_region'),#all available listFilters(ensembl), search for filters searchFilters(mart=ensembl_db,pattern='sequence')
                 values=list(genomic_region),#all available listAttributes(ensembl), search for attributes searchAtrributes(mart=ensembl_db,pattern='sequence')
                 mart=ensembl_db)
  return(ensembl)
}

genes<-sapply(genomic_regions,ensembl_query)
genes<-unlist(genes)

#create dataframe from stacked named lists
df<-stack(genes)

#write output csv
write.csv(df,sprintf('/Volumes/motilin/balint/numt/mice_strains_numt/gene_ids/%s.csv',current_strain))

