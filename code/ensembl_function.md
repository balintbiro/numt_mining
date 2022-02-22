intragenic_numts
================
Balint Biro
2/22/2022

``` r
#documentation -->https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
#load the package
library(biomaRt)
```

``` r
#create vector from mice strains
mice_strains <- c('mmusculus','mmwsbeij','mmbalbcj','mmcasteij','mmpwkphj','mmfvbnj','mmdba2j',
                  'mm129s1svimj','mmc3hhej','mmnodshiltj','mmnzohlltj','mmc57bl6nj',
                  'mmcbaj','mmaj','mmakrj','mmlpj','mspretus')
```

Create function for querying Ensembl. First step, select the desired
database, all available listEnsembl(). Select the desired dataset, all
available listDatasets(ensembl). Select the mirror, options
–>useast,uswest,asia,www

``` r
ensembl_query <- function(mice_strain){
  genomic_regions<-scan(
    sprintf
    (
      '/Volumes/motilin/balint/numt/mice_strains_numt/genomic_regions/%s_gene_ensembl.txt',
            mice_strain
    ),character())
  not_strains <- c('mmusculus','mspretus')
  if(mice_strain %in% not_strains){
    ensembl_db<-useEnsembl(biomart='genes',
                         dataset=sprintf('%s_gene_ensembl',mice_strain),
                         mirror='uswest')
  }else{
    ensembl_db<-useEnsembl(biomart='mouse_strains',
                         dataset=sprintf('%s_gene_ensembl',mice_strain),
                         mirror='uswest')
  }
  ensembl<-getBM(attributes=c('external_gene_name','goslim_goa_description','name_1006'),
                 filters=c('chromosomal_region'),
                 values=genomic_regions,
                 mart=ensembl_db)
  return(ensembl)
}
gos <- sapply(mice_strains,ensembl_query)
```

``` r
mice_strains
```

    ##  [1] "mmusculus"    "mmwsbeij"     "mmbalbcj"     "mmcasteij"    "mmpwkphj"    
    ##  [6] "mmfvbnj"      "mmdba2j"      "mm129s1svimj" "mmc3hhej"     "mmnodshiltj" 
    ## [11] "mmnzohlltj"   "mmc57bl6nj"   "mmcbaj"       "mmaj"         "mmakrj"      
    ## [16] "mmlpj"        "mspretus"

``` r
actual_strain <- mice_strains[1]
mmusculus <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[2]
mmwsbeij <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[3]
mmbalbcj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[4]
mmcasteij <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[5]
mmpwkphj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[6]
mmfvbnj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[7]
mmdba2j <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[8]
mm129s1svimj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[9]
mmc3hhej <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[10]
mmnodshiltj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[11]
mmnzohlltj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[12]
mmc57bl6nj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[13]
mmcbaj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[14]
mmaj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[15]
mmakrj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[16]
mmlpj <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))

actual_strain <- mice_strains[17]
mspretus <- as.data.frame(cbind(gene_name=unlist(gos['external_gene_name',actual_strain]),
                      go_description=unlist(gos['goslim_goa_description',actual_strain]),
                      go_name=unlist(gos['name_1006',actual_strain])))
```

``` r
dfs <- list(mmusculus,mmwsbeij,mmbalbcj,mmcasteij,mmpwkphj,mmfvbnj,mmdba2j,
                  mm129s1svimj,mmc3hhej,mmnodshiltj,mmnzohlltj,mmc57bl6nj,
                  mmcbaj,mmaj,mmakrj,mmlpj,mspretus)
names(dfs) <- c('mmusculus','mmwsbeij','mmbalbcj','mmcasteij','mmpwkphj','mmfvbnj','mmdba2j',
                  'mm129s1svimj','mmc3hhej','mmnodshiltj','mmnzohlltj','mmc57bl6nj',
                  'mmcbaj','mmaj','mmakrj','mmlpj','mspretus')
```

Function for writing outputs.

``` r
write_output <- function(df){
  tracker <<- tracker+1
  write.csv(dfs[tracker],
            sprintf('/Volumes/motilin/balint/numt/trials/%s.csv',
                   mice_strains[tracker]))
}
tracker <- 0
sapply(dfs,
       write_output)
```

    ## $mmusculus
    ## NULL
    ## 
    ## $mmwsbeij
    ## NULL
    ## 
    ## $mmbalbcj
    ## NULL
    ## 
    ## $mmcasteij
    ## NULL
    ## 
    ## $mmpwkphj
    ## NULL
    ## 
    ## $mmfvbnj
    ## NULL
    ## 
    ## $mmdba2j
    ## NULL
    ## 
    ## $mm129s1svimj
    ## NULL
    ## 
    ## $mmc3hhej
    ## NULL
    ## 
    ## $mmnodshiltj
    ## NULL
    ## 
    ## $mmnzohlltj
    ## NULL
    ## 
    ## $mmc57bl6nj
    ## NULL
    ## 
    ## $mmcbaj
    ## NULL
    ## 
    ## $mmaj
    ## NULL
    ## 
    ## $mmakrj
    ## NULL
    ## 
    ## $mmlpj
    ## NULL
    ## 
    ## $mspretus
    ## NULL
