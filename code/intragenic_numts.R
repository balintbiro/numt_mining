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
      values=list(chromosome='Mt',start='3000',end='8800'),
      mart=ensembl)

data<-read.csv('Z:/balint/numt/mice_strains_numts.csv')
subdf<-head(data)
proba<-function(row){
  return(c(row['g_length'],row['mt_length']))
}
ensembl_query<-function(row){
  strain<-strsplit(as.character(row['organism_name']),'_',FALSE)[[1]][3]
  chr_name<-row['g_id']
  numt_start<-as.numeric(row['g_start'])
  numt_end<-numt_start+as.numeric(row['g_length'])
  
  ensembl<-useEnsembl(biomart='mouse_strains',#select the desired database
                      dataset=paste('mm',strain,'_gene_ensembl',sep=''),#select the desired dataset
                      mirror='asia')#select the mirror, options -->useast,uswest,asia,www
  informations<-getBM(attributes=c('external_gene_name','uniprotswissprot'),
               filters=c('chromosome_name','start','end'),
               values=list(chromosome=chr_name,start=numt_start,end=numt_end),
               mart=ensembl)
  return(informations)
}

para<-apply(data[1,],1,ensembl_query)
