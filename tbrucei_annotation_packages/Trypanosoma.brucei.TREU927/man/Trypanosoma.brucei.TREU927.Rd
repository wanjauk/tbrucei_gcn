\name{Trypanosoma.brucei.TREU927}
\docType{data}
\alias{Trypanosoma.brucei.TREU927}
\title{Annotation package that collates several annotation resources.}
\description{This data object was automatically created by Keith Hughitt <khughitt@umd.edu>. It
  represents a collection of annotation packages that can be used as a
  single object named Trypanosoma.brucei.TREU927.  This object can be used with the
  standard four accessor method for all AnnotationDbi objects.  Namely:
  \code{cols}, \code{keytype}, \code{keys} and \code{select}.  Users are
  encouraged to read the vignette from the \code{OrganismDbi} package for
  more details.}
\usage{library(Trypanosoma.brucei.TREU927)}
\keyword{datasets}
\examples{
  Trypanosoma.brucei.TREU927
  cls <- columns(Trypanosoma.brucei.TREU927)
  cls
  cls <- cls[c(1,19,45)]
  kts <- keytypes(Trypanosoma.brucei.TREU927)
  kt <- kts[2]
  kts
  ks <- head(keys(Trypanosoma.brucei.TREU927, keytype=kts[2]))
  ks
  res <- select(Trypanosoma.brucei.TREU927, keys=ks, columns=cls, keytype=kt)
  head(res)
}
