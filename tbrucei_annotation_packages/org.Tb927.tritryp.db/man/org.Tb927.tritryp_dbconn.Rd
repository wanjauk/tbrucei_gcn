
\name{org.Tb927.tritryp_dbconn}

\alias{org.Tb927.tritryp_dbconn}
\alias{org.Tb927.tritryp_dbfile}
\alias{org.Tb927.tritryp_dbschema}
\alias{org.Tb927.tritryp_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  org.Tb927.tritryp_dbconn()
  org.Tb927.tritryp_dbfile()
  org.Tb927.tritryp_dbschema(file="", show.indices=FALSE)
  org.Tb927.tritryp_dbInfo()
}

\arguments{
  \item{file}{
    A connection, or a character string naming the file to print to (see
    the \code{file} argument of the \code{\link[base]{cat}} function for
    the details).
  }
  \item{show.indices}{
    The CREATE INDEX statements are not shown by default.
    Use \code{show.indices=TRUE} to get them.
  }
}

\details{
  \code{org.Tb927.tritryp_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{org.Tb927.tritryp_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{org.Tb927.tritryp_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{org.Tb927.tritryp_dbschema} prints the schema definition of the
  package annotation DB.

  \code{org.Tb927.tritryp_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{org.Tb927.tritryp_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{org.Tb927.tritryp_dbfile}: a character string with the path to the
  package annotation DB.

  \code{org.Tb927.tritryp_dbschema}: none (invisible \code{NULL}).

  \code{org.Tb927.tritryp_dbInfo}: none (invisible \code{NULL}).
}

\seealso{
  \code{\link[DBI:dbSendQuery]{dbGetQuery}},
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbconn}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbfile}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbschema}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbInfo}}
}

\examples{
  library(DBI)
  ## Count the number of rows in the "genes" table:
  dbGetQuery(org.Tb927.tritryp_dbconn(), "SELECT COUNT(*) FROM genes")

  org.Tb927.tritryp_dbschema()

  org.Tb927.tritryp_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

