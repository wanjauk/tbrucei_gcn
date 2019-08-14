# TxDb.TbruceiTREU927.tritryp43.genes

Transcript annotation package for *Trypanosoma brucei brucei TREU927*, based on
annotated genes from [TriTrypDB 43](http://tritrypdb.org/tritrypdb/).

This package was generated using the tools from
[https://github.com/elsayed-lab/eupathdb-organismdb](github.com/eupathdb-organismdb).

Installation
------------

You can install the latest version from Github using:

``` r
library('devtools')
install_github('elsayed-lab/TxDb.TbruceiTREU927.tritryp43.genes')
```

Usage
-----

This package is based on the Bioconductor
[AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
interface. As such, the methods for interacting with this package are similar
to the ways one can interact with other commonly-used annotation packages such as
[TxDb.Hsapiens.UCSC.hg19.knownGene](http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html).

Example usage:

```r
library(TxDb.TbruceiTREU927.tritryp43.genes)

# list available fields to query
columns(TxDb.TbruceiTREU927.tritryp43.genes)

# get first 10 genes
gene_ids = head(keys(TxDb.TbruceiTREU927.tritryp43.genes), 10)

# gene coordinates and strand
genes = AnnotationDbi::select(TxDb.TbruceiTREU927.tritryp43.genes, 
                              keys=gene_ids, 
                              keytype='GENEID', 
                              columns=c('TXSTART', 'TXEND', 'TXSTRAND'))

head(genes)
```

For more information, check out the [AnnotationDbi - Introduction to Annotation
packages vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

Additional resources that may be helpful:

1. http://www.bioconductor.org/help/workflows/annotation-data/
2. http://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
3. http://training.bioinformatics.ucdavis.edu/docs/2012/05/DAV/lectures/annotation/annotation.html
