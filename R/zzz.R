.onAttach <- function(libname, pkgname) {
    v <- read.dcf(system.file("DESCRIPTION", package = pkgname), fields = "Version")
    msg <- paste0("Loading ", pkgname, " v", v[[1]])
    packageStartupMessage(msg)
}

utils::globalVariables(c(
    "Dim1", "Dim2", "passed", "fdr", "threshold"
))