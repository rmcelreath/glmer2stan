.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  theLib <- dirname(system.file(package = "glmer2stan"))
  pkgdesc <- packageDescription("glmer2stan", lib.loc = theLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  packageStartupMessage(paste("glmer2stan (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
} 