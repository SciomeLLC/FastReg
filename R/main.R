#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param phenotype bin.resp | num.resp. Default: bin.resp
#' @param regression.type logistic | linear. Default: logistic
#' @param Pvalue.dist t.dist | norm.dist. Default: t.dist
#' @param output.exclude.covar TRUE | FALSE. Default: TRUE. Excludes covariates in the results output
#' @param maf.threshold Default: 1e-13. Minor Allele frequency threshold. Acceptable values: 0-0.5
#' @param hwe.threshold Default: 1e-13. Hardy-Weinberg Equilibrium threshold. Acceptable values: 0-0.5
#' @param no.intercept TRUE | FALSE. Default: FALSE.
#' @param colinearity.rsq Default: 1.0. Filter covariates if above R-squared colinearity threshold. Acceptage values: 0.8-1.0
#' @param poi.block.size Default: 0. Overrides auto-calculated poi block size.
#' @param max.iter Default: 6. Number of logistic regression iterations. Doesn't apply when regression.type = linear.
#' @param rel.conv.tolerance Default: 0.01. Relative convergence threshold for POIs.
#' @param abs.conv.tolerance Default: 0.01. Absolute convergence threshold for POIs.
#' @param max.threads Default: 0. Overrides number of threads used. By default FastReg will use the max number of threads available - 1.
#' @param pheno.file relative path to phenotype file. Contains outcomes that will be modelled.
#' @param pheno.rowname.cols Default: ind. column to be treated as the subject identifier when matching across files. Can be multiple columns (comma separated).
#' @param pheno.file.delim tab|space|comma. Default: tab. delimiter used in phenotype file.
#' @param covar.file relative path to covariates file.
#' @param covar.rowname.cols Default: ind. column to be treated as the subject identifier when matching across files. Can be multiple columns (comma separated).
#' @param covar.file.delim tab|space|comma. Default: tab. delimiter used in covariate file.
#' @param POI.file.dir relative path to the directory with POI file(s).
#' @param POI.file.delim tab|space|comma. Default: tab. delimiter used in POI file if POI.file.format != h5
#' @param POI.file.format h5|txt. Default: h5. POI data file format.
#' @param POI.type genotype|dosage. Default: dosage.
#' @param POI.effect.type additive|dominant|recessive. Default: additive.
#' @param covariates Default: c("age","sex","eth","treatment","severity"). List of column names of the covariates file.
#' @param covariate.type Default: c("numeric","categorical","categorical","categorical","categorical"). List of column data type for each covariate listed in covariates.
#' @param covariate.standardize Default: c(TRUE,FALSE,FALSE,FALSE,FALSE). List of booleans representing which covariates are to be standardized.
#' @param covariate.levels Default: c("", "F,M","afr,asi,eur","Placebo,Test","Very Low,Low,Moderate,High,Very High,Extreme"). List of unique values for categorical columns. (comma separated for each unique value)
#' @param covariate.ref.level Default: c("","F","eur","Placebo","Very Low"). List of reference level for covariate value (comma separated)
#' @param POI.covar.interactions Default: c(""). List of poi-covariate interactions
#' @param split.by Default: c(""). List of covariates values to split by when stratifying the data.
#' @param output.dir Default: results. Relative path to directory for output. Directory will be created if it doesn't exist
#' @param compress.results TRUE|FALSE. Default: FALSE. Compress results inside the output directory.
#' @param config.file deprecated. No effect.
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import Rcpp
#' @import RcppArmadillo
#' @import stats
#' @import parallel
#' @import data.table

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
    POI.file.dir = "testdata/",
    POI.file.delim = "tab",
    POI.file.format = "h5",
    POI.type = "dosage",
    POI.effect.type = "additive",
    covariates = c("age", "sex", "eth", "treatment", "severity"),
    covariate.type = c("numeric", "categorical", "categorical", "categorical", "categorical"),
    covariate.standardize = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    covariate.levels = c("", "F,M", "afr,asi,eur", "Placebo,Test", "Very Low,Low,Moderate,High,Very High,Extreme"),
    covariate.ref.level = c("", "F", "eur", "Placebo", "Very Low"),
    POI.covar.interactions = c(""),
    split.by = c(""),
    output.dir = "test",
    compress.results = FALSE) {
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
    POI.file.dir,
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

#' hdf5convert a function to convert textual data to hdf5 format supported by FastReg()
#' @param dataFile file name of the text data
#' @param h5Dir file dir to write hdf5 files to
#' @param headerRow row containing the header
#' @param idCol column with the ids
#' @param dataCol start of the data columns in the text data
#' @param buffSize buffer size
#' @param transpose boolean to transpose the data
#' @param chunkEdge chunk size for the hdf5 file
#' @param vcf TRUE or FALSE
#' @param delimiter delimiter used in the dataFile
#' @param gz bool to gzip the h5 dataset or not
#' @param poiPerFile =-1
#' @param singleFile =FALSE
#' @param serverThreads =-1
#' @param serverMem =-1.0
#' @return Boolean stating if the result file was created succesfully
#' @export

fastR_hdf5convert <- function(
    dataFile,
    h5Dir,
    headerRow = 1,
    idCol = 1,
    dataCol = 2,
    buffSize = 1.0,
    transpose = FALSE,
    chunkEdge = 100,
    vcf = FALSE,
    delimiter = " \t",
    gz = FALSE,
    poiPerFile = -1,
    singleFile = FALSE,
    serverThreads = -1,
    serverMem = -1.0) {
  if (!is.character(dataFile) || length(dataFile) > 1) {
    cat("Error: dataFile must be a single-member character vector\n")
    return(FALSE)
  }
  if (!is.character(h5Dir) || length(h5Dir) > 1) {
    cat("Error: h5Dir must be a single-member character vector\n")
    return(FALSE)
  } else if (file.exists(h5Dir)) {
    cat("Error: file or folder named \"", h5Dir, "\" already exists\n", sep = "")
    return(FALSE)
  }
  dir.create(h5Dir)
  if (!is.numeric(headerRow) || headerRow != as.integer(headerRow) || length(headerRow) > 1 || headerRow < 1) {
    cat("Error: headerRow must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  headerRow <- as.integer(headerRow)
  if (!is.numeric(idCol) || idCol != as.integer(idCol) || length(idCol) > 1 || idCol < 1) {
    cat("Error: idCol must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  idCol <- as.integer(idCol)
  if (!is.numeric(dataCol) || dataCol != as.integer(dataCol) || length(dataCol) > 1 || dataCol <= idCol) {
    cat("Error: headerRow must be a single-member integer vector with value > idCol\n")
    return(FALSE)
  }
  dataCol <- as.integer(dataCol)
  if (!is.numeric(buffSize) || length(buffSize) > 1 || buffSize <= 0) {
    cat("Error: buffSize must be a single-member numeric vector with value > 0\n")
    return(FALSE)
  }
  buffSize <- as.double(buffSize)
  if (!is.logical(transpose) || length(transpose) > 1) {
    cat("Error: transpose must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(chunkEdge) || chunkEdge != as.integer(chunkEdge) || length(chunkEdge) > 1 || chunkEdge < 1) {
    cat("Error: chunkEdge must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  chunkEdge <- as.integer(chunkEdge)
  if (!is.logical(vcf) || length(vcf) > 1) {
    cat("Error: vcf must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.character(delimiter) || length(delimiter) > 1) {
    cat("Error: delimiter must be a single-member character vector\n")
    return(FALSE)
  }
  if (!is.logical(gz) || length(gz) > 1) {
    cat("Error: gz must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(poiPerFile) || poiPerFile != as.integer(poiPerFile) || length(poiPerFile) > 1 || poiPerFile < chunkEdge) {
    if (poiPerFile != -1) {
      cat("Error: poiPerFile must be a single-member integer vector with value >= chunkEdge\n")
      return(FALSE)
    }
  }
  poiPerFile <- as.integer(poiPerFile)
  if (!is.logical(singleFile) || length(singleFile) > 1) {
    cat("Error: singleFile must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(serverThreads) || serverThreads != as.integer(serverThreads) || length(serverThreads) > 1) {
    cat("Error: serverThreads must be a single-member integer vector\n")
    return(FALSE)
  }
  serverThreads <- as.integer(serverThreads)
  if (!is.numeric(serverMem) || length(serverMem) > 1) {
    cat("Error: serverMem must be a single-member numeric vector\n")
    return(FALSE)
  }
  serverMem <- as.double(serverMem)
  fastR_hdf5convert_cpp(
    dataFile,
    h5Dir,
    headerRow,
    idCol,
    dataCol,
    buffSize,
    transpose,
    chunkEdge,
    vcf,
    delimiter,
    gz,
    poiPerFile,
    singleFile,
    serverThreads,
    serverMem
  )
}