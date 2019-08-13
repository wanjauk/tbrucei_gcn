\name{TxDb.TbruceiTREU927.tritryp43.genes}
\docType{package}

\alias{TxDb.TbruceiTREU927.tritryp43.genes-package}
\alias{TxDb.TbruceiTREU927.tritryp43.genes}
\alias{TxDb.TbruceiTREU927.tritryp43.genes}


\title{Annotation package for TxDb object(s)}

\description{  
  This package loads one or more TxDb objects.  Such TxDb
  objects are an R interface to prefabricated databases contained by
  this package.

  The names of any objects exposed by this package indicate the origin and
  resources exposed.  So for example TxDb.Hsapiens.UCSC.hg19.knownGene
  would be a TxDb object, of Homo sapiens data from UCSC build
  hg19 based on the knownGene Track.
}

\note{
  This data package was made from resources at TriTrypDB on
  2019-08-06 14:43:34 +0300 (Tue, 06 Aug 2019) and based on the 43
}

\author{Keith Hughitt <khughitt@umd.edu> [cre]}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(TxDb.TbruceiTREU927.tritryp43.genes)
## list the contents that are loaded into memory
ls('package:TxDb.TbruceiTREU927.tritryp43.genes')
## show the db object that is loaded by calling it's name
TxDb.TbruceiTREU927.tritryp43.genes

}

\keyword{package}
\keyword{data}
