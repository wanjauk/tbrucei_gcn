datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Tb927.tritryp <- function() showQCData("org.Tb927.tritryp", datacache)
org.Tb927.tritryp_dbconn <- function() dbconn(datacache)
org.Tb927.tritryp_dbfile <- function() dbfile(datacache)
org.Tb927.tritryp_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Tb927.tritryp_dbInfo <- function() dbInfo(datacache)

org.Tb927.tritrypORGANISM <- "Trypanosoma brucei"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Tb927.tritryp.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Tb927.tritryp.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Tb927.tritryp_dbconn())
}

