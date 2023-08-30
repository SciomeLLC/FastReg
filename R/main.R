#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import Rcpp
#' @import RcppArmadillo

FastReg <- function(
  phenotype = "bin.resp",
  regression.type = "logistic",
  Pvalue.dist = "t.dist",
  output.exclude.covar = TRUE,
  maf.threshold = 1e-13,
  hwe.threshold = 1e-13,
  no.intercept = FALSE,
  colinearity.rsq = 1.0,
  poi.block.size = 0,
  max.iter = 6,
  rel.conv.tolerance = 0.01,
  abs.conv.tolerance = 0.01,
  max.threads = 0,
  pheno.file = "testdata_1k_by_5k.bin.pheno.txt",
  pheno.rowname.cols = "ind",
  pheno.file.delim = "tab",
  covar.file = "testdata_1k_by_5k.covar.txt",
  covar.rowname.cols = "ind",
  covar.file.delim = "tab",
  POI.file = "testdata_1k_by_5k.poi.h5",
  POI.file.delim = "tab",
  POI.file.format = "h5",
  POI.type = "dosage",
  POI.effect.type = "additive",
  covariates = c("age","sex","eth","treatment","severity"),
  covariate.type = c("numeric","categorical","categorical","categorical","categorical"),
  covariate.standardize = c(TRUE,FALSE,FALSE,FALSE,FALSE),
  covariate.levels = c("", "F,M","afr,asi,eur","Placebo,Test","Very Low,Low,Moderate,High,Very High,Extreme"),
  covariate.ref.level = c("","F","eur","Placebo","Very Low"),
  POI.covar.interactions = c(""),
  split.by = c(""),
  output.dir = "test",
  compress.results = FALSE
  ) {
  FastRegCpp(
    phenotype,
    regression.type,
    Pvalue.dist,
    output.exclude.covar,
    maf.threshold,
    hwe.threshold,
    no.intercept,
    colinearity.rsq,
    poi.block.size,
    max.iter,
    rel.conv.tolerance,
    abs.conv.tolerance,
    max.threads,
    pheno.file,
    pheno.rowname.cols,
    pheno.file.delim,
    covar.file,
    covar.rowname.cols,
    covar.file.delim,
    POI.file,
    POI.file.delim,
    POI.file.format,
    POI.type,
    POI.effect.type,
    covariates,
    covariate.type,
    covariate.standardize,
    covariate.levels,
    covariate.ref.level,
    POI.covar.interactions,
    split.by,
    output.dir,
    compress.results
  )
}
