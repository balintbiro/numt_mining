
#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#load the package
library(biomaRt)

#define function for creating ENsembl database and dataset
create_ensembl_db_ds<-function(mouse_strain){
  strain<-strsplit(mouse_strain,'_',FALSE)
  ensembl_strain<-paste('mm',strain[[1]][3],'_gene_ensembl',sep='')
  ensembl<-useEnsembl(biomart='mouse_strains',#select the desired database, all available listEnsembl()
                      dataset=ensembl_strain,#select the desired dataset, all available listDatasets(ensembl)
                      mirror='asia')#select the mirror, options -->useast,uswest,asia,www
  return(ensembl)
}

ensembl_db<-sapply(unique(data['organism_name'][[1]]),create_ensembl_db_ds)

#define function for querying Ensembl
ensembl_query<-function(row){
  strain<-strsplit(as.character(row['organism_name']),'_',FALSE)[[1]][3]
  chr_name<-row['g_id']
  numt_start<-as.numeric(row['g_start'])
  numt_end<-numt_start+as.numeric(row['g_length'])
  
  #check the avaialable databases
  #listEnsembl()
  #check the avaiable datasets
  #listDatasets(ensembl)
  ensembl<-useEnsembl(biomart='mouse_strains',#select the desired database
                      dataset=paste('mm',strain,'_gene_ensembl',sep=''),#select the desired dataset
                      mirror='asia')#select the mirror, options -->useast,uswest,asia,www
  
  #check the available filters and attributes
  #listFilters(ensembl)
  #listAttributes(ensembl)
  informations<-getBM(attributes=c('external_gene_name','uniprotswissprot'),
               filters=c('chromosome_name','start','end'),
               values=list(chromosome=chr_name,start=numt_start,end=numt_end),
               mart=ensembl)
  return(as.vector(unname(informations)))
}

#load the csv file
data<-read.csv('Z:/balint/numt/ensembl_tester.csv')

para<-apply(data,1,ensembl_query)
write.csv(para,'Z:/balint/numt/proba_output.csv')
ensembl <- useEnsembl(biomart='mouse_strains',#select the desired database
                      dataset='mm129s1svimj_gene_ensembl',#select the desired dataset
                      mirror='asia')#select the mirror, options -->useast,uswest,asia,www







































