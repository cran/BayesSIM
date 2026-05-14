# R/zzz.R
#' @keywords internal
"_PACKAGE"
.onLoad <- function(lib, pkg) {
  loadNamespace("nimble")
}
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading BayesSIM...\n")
    }


