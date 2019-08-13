# org.Tb927.tritryp.db

Genome-wide annotation package for *Trypanosoma brucei brucei TREU927*, based on
annotations from [TriTrypDB 43](http://tritrypdb.org/tritrypdb/).

This package was generated using the tools from
[https://github.com/elsayed-lab/eupathdb-organismdb](github.com/eupathdb-organismdb).

Installation
------------

You can install the latest version from Github using:

``` r
library('devtools')
install_github('elsayed-lab/org.Tb927.tritryp.db')
```

Usage
-----

This package is based on the Bioconductor
[AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
interface. As such, the methods for interacting with this package are similar
to the ways one can interact with other commonly-used annotation packages such as
[org.Hs.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html).

Example usage:

```r
library(org.Tb927.tritryp.db)

# list available fields to query
columns(org.Tb927.tritryp.db)

# get first 10 genes
gene_ids = head(keys(org.Tb927.tritryp.db), 10)

# gene names and descriptions
annotations = AnnotationDbi::select(org.Tb927.tritryp.db, 
                                    keys=gene_ids, 
                                    keytype='GID', 
                                    columns=c('CHROMOSOME', 'GENEDESCRIPTION'))
head(annotations)

# GO terms
go_terms = AnnotationDbi::select(org.Tb927.tritryp.db, 
                                 keys=gene_ids, 
                                 keytype='GID', 
                                 columns=c('GO', 'ONTOLOGYALL'))
head(go_terms)

# KEGG pathways
kegg_paths = AnnotationDbi::select(org.Tb927.tritryp.db,
                                   keys=gene_ids, 
                                   keytype='GID', 
                                   columns=c('KEGG_NAME', 'KEGG_PATH'))
head(kegg_paths)
```

For more information, check out the [AnnotationDbi - Introduction to Annotation
packages vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf).

Additional resources that may be helpful:

1. http://www.bioconductor.org/help/workflows/annotation-data/
2. http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
3. http://training.bioinformatics.ucdavis.edu/docs/2012/05/DAV/lectures/annotation/annotation.html
