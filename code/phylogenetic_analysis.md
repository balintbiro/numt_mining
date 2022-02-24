phylogenetics
================
Balint Biro
2/23/2022

Documnetations: - <https://yulab-smu.top/treedata-book/chapter4.html> -
<https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html>
- <https://rdrr.io/cran/phangorn/f/vignettes/Trees.Rmd> Related
articles: - <https://academic.oup.com/jhered/article/98/3/243/2188032> -
<https://academic.oup.com/sysbio/article/54/6/952/1630741>

``` r
#loading the required packages
library(seqinr)
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following objects are masked from 'package:seqinr':
    ## 
    ##     as.alignment, consensus

``` r
library(adegenet)
```

    ## Loading required package: ade4

    ## 
    ##    /// adegenet 2.1.5 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

    ## 
    ## Attaching package: 'adegenet'

    ## The following object is masked from 'package:phangorn':
    ## 
    ##     AICc

``` r
library(ape)
library(ggtree)
```

    ## ggtree v3.2.1  For help: https://yulab-smu.top/treedata-book/
    ## 
    ## If you use ggtree in published research, please cite the most appropriate paper(s):
    ## 
    ## 1. Guangchuang Yu. Using ggtree to visualize data on tree-like structures. Current Protocols in Bioinformatics. 2020, 69:e96. doi:10.1002/cpbi.96
    ## 2. Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. Molecular Biology and Evolution. 2018, 35(12):3041-3043. doi:10.1093/molbev/msy194
    ## 3. Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628

    ## 
    ## Attaching package: 'ggtree'

    ## The following object is masked from 'package:ape':
    ## 
    ##     rotate

``` r
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:ade4':
    ## 
    ##     score

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:ggtree':
    ## 
    ##     expand

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:ggtree':
    ## 
    ##     collapse

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:ape':
    ## 
    ##     complement

    ## The following object is masked from 'package:seqinr':
    ## 
    ##     translate

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(ggplot2)
```

UPGMA tree with bootstrap support

``` r
pseudogenes <- read.phyDat('../data/pseudogenes.dna',
                          format='interleaved')

set.seed(1)
fun <- function(x) upgma(dist.ml(x,model='JC69'))
bs_upgma <- bootstrap.phyDat(pseudogenes,  fun, bs=100)
#png('../results/bs_upgma.png',width=800, height=600)
ggtree(root(consensus(bs_upgma),outgroup='rat_CYTB'))+
  geom_tiplab(size=5,align=TRUE)+
  xlim(0,5.5)+
  geom_cladelabel(node=37, label="numt", 
                  color="dark green",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.5)+
  geom_cladelabel(node=38, label="CYTB", 
                  color="brown",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.5)+
  geom_cladelabel(node=10, label="outgroup", 
                  color="blue",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.5)
```

![](Untitledphylogenetic_analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#dev.off()
```

Parsimony tree with bootstrap support

``` r
pseudogenes <- read.phyDat('../data/pseudogenes.dna',
                          format='interleaved')
treepars <- pratchet(pseudogenes,
                     trace=0)
#png('../results/bs_mp.png',width=800, height=600)
ggtree(root(treepars,outgroup='rat_CYTB'))+
  geom_tiplab(align=TRUE,size=5)+
  xlim(0,7.8)+
  geom_cladelabel(node=38, label="numt", 
                  color="dark green",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.85)+
  geom_cladelabel(node=39, label="CYTB", 
                  color="brown",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.85)+
  geom_cladelabel(node=4, label="outgroup", 
                  color="blue",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=0.85)
```

![](Untitledphylogenetic_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#dev.off()
```

Maximum Likelihood tree

``` r
pseudogenes <- read.phyDat('../data/pseudogenes.dna',
                          format='interleaved')
dm <- dist.ml(pseudogenes,
              model='JC69')
treeupgma <- root(upgma(dm),outgroup='rat_CYTB')
fit <- pml(treeupgma,data=pseudogenes)
#png('../results/ml.png',width=800, height=600)
ggtree(root(fit$tree,outgroup='rat_CYTB'))+
  geom_tiplab(aling=TRUE,size=5)+
  xlim(0,0.3)+
  geom_cladelabel(node=37, label="numt", 
                  color="dark green",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE,offset=.035)+
  geom_cladelabel(node=38, label="CYTB", 
                  color="brown",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset=-.115)+
  geom_cladelabel(node=18, label="outgroup", 
                  color="blue",
                  fontsize=5, alpha=0.5,
                  barsize=2, hjust=-.1,
                  align=TRUE, offset.text = -.12)
```

![](Untitledphylogenetic_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#dev.off()
```

``` r
#bs consensus upgma tree
pseudogenes <- read.phyDat('../data/pseudogenes.dna',
                          format='interleaved')

dm <- dist.ml(pseudogenes,
        model='JC69')

fun <- function(x) upgma(dist.ml(x,model='JC69'))
bs_upgma <- bootstrap.phyDat(pseudogenes,  fun, bs=100)

treeupgma <- upgma(dm)
#png('../results/bs_trees.png',width=800, height=600)
ggdensitree(bs_upgma, alpha=0.15, color='lightblue',align.tips=TRUE, size=1.5)+
  geom_tiplab(size=5, offset=0.003)+
  xlim(-0.2,0.05)+
  geom_cladelabel(node=37, label="CYTB", 
                  color="brown",
                  fontsize=8, alpha=0.5,
                  barsize=2, offset=0.03, offset.text=0.005,
                  angle=90)+
  geom_cladelabel(node=36, label="numt", 
                  color="darkgreen",
                  fontsize=8, alpha=0.5,
                  barsize=2, offset=0.03, offset.text=0.005,
                  angle=90)
```

    ## Warning in max(sp.df$x, na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in min(y): no non-missing arguments to min; returning Inf

    ## Warning in max(y): no non-missing arguments to max; returning -Inf

    ## Warning in max(sp.df$x, na.rm = TRUE): no non-missing arguments to max;
    ## returning -Inf

    ## Warning in min(y): no non-missing arguments to min; returning Inf

    ## Warning in max(y): no non-missing arguments to max; returning -Inf

![](Untitledphylogenetic_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#dev.off()
```
