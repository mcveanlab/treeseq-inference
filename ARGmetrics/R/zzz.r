## for NAMESPACE exporting
#' @export
kc.dist <- NULL

.onLoad <- function(lib, pkg, ...) {
  if (suppressWarnings(require("treespace", character.only = TRUE))) {
    kc.dist <<- treespace::treeDist
  } else {
    warning("Package 'treespace' not installed, using built-in alternative to treeDist()")
    kc.dist <<- treeDist
  }
}

#' @useDynLib ARGmetrics
#' @importFrom Rcpp sourceCpp
NULL
