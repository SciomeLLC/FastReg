#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import Rcpp
#' @import RcppArmadillo

FastReg <- function(config.file, ...) {
  if(is.null(config.file)) stop("Configuration file path must be specified");
  FastRegCpp(config.file)
}
