#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import Rcpp
#' @import RcppArmadillo

FastReg <- function(
  # phenotype = "bin.resp",
  # regression.type = "logistic",
  # Pvalue.dist = "t.dist",
  # output.exclude.covar = TRUE,
  # maf.threshold = 1e-13,
  # hwe.threshold = 1e-13,
  # no.intercept = FALSE,
  # colinearity.rsq = 1.0,
  # poi.block.size = 0,
  # max.iter = 6,
  # rel.conv.tolerance = 0.01,
  # abs.conv.tolerance = 0.01,
  # max.threads = 0,
  
  config.file
  ) {
  if(is.null(config.file)) stop("Configuration file path must be specified");
  FastRegCpp(config.file)
}
