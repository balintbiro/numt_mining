
#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#load the package
library(biomaRt)

#load the csv file
data<-read.csv('Z:/balint/numt/ensembl_tester.csv')

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
  
  mouse_strain<-row['organism_name']
  ensembl<-ensembl_db[mouse_strain][[1]]
  chr_name<-row['g_id']
  numt_start<-as.numeric(row['g_start'])
  numt_end<-numt_start+as.numeric(row['g_length'])
  
  informations<-getBM(attributes=c('external_gene_name','uniprotswissprot'),
               filters=c('chromosome_name','start','end'),#all available listFilters(ensembl)
               values=list(chromosome=chr_name,start=numt_start,end=numt_end),#all available listAttributes(ensembl)
               mart=ensembl)
  if (is.logical(informations[[1]])){
    gene_names<<-append(gene_names,NA)
    uniprot_ids<<-append(uniprot_ids,NA)
  }
  else{
    gene_names<<-append(gene_names,list(informations['external_gene_name'][[1]]))
    uniprot_ids<<-append(uniprot_ids,list(informations['uniprotswissprot'][[1]]))
  }
}

#define function for querying Ensembl
ensembl_query<-function(row){
  out<-tryCatch(
    
    #try branch
    {
      mouse_strain<-row['organism_name']
      ensembl<-ensembl_db[mouse_strain][[1]]
      chr_name<-row['g_id']
      numt_start<-as.numeric(row['g_start'])
      numt_end<-numt_start+as.numeric(row['g_length'])
      
      informations<-getBM(attributes=c('external_gene_name','uniprotswissprot'),
                          filters=c('chromosome_name','start','end'),#all available listFilters(ensembl)
                          values=list(chromosome=chr_name,start=numt_start,end=numt_end),#all available listAttributes(ensembl)
                          mart=ensembl)
      if (is.logical(informations[[1]])){
        gene_names<<-append(gene_names,NA)
        uniprot_ids<<-append(uniprot_ids,NA)
      }
      else{
        gene_names<<-append(gene_names,list(informations['external_gene_name'][[1]]))
        uniprot_ids<<-append(uniprot_ids,list(informations['uniprotswissprot'][[1]]))
      }
    },
    
    #except branch
    error=function(row){
      gene_names<<-append(gene_names,NA)
      uniprot_ids<<-append(uniprot_ids,NA)
    }
  )
  
}

gene_names<-c()
uniprot_ids<-c()
proba<-apply(data,1,ensembl_query)


