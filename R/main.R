#' FastReg a function to perform fast simple linear or logistic regression
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
#' @param max.openmp.threads Default: 2. Overrides number of threads used. By default FastReg will use 2 openmp threads on Ubuntu otherwise 1.
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
#' @param max.workers Default: 0. Number of processes FastReg will spawn. By default FastReg will spawn min(# of poi files, # of available workers)
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import Rcpp
#' @import RcppArmadillo
#' @import parallel
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
    max.openmp.threads = 2,
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
    compress.results = FALSE,
    max.workers = 0) {
  if (max.openmp.threads <= 0) {
    cat("Error: max.openmp.threads must be a positive integer.\n")
    return(FALSE)
  }

  if (max.workers < 0) {
    cat("Error: max.workers must be a positive integer.\n")
    return(FALSE)
  }
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
    max.openmp.threads,
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
    compress.results,
    max.workers
  )
}

#' FastRegImport a function to convert textual data to hdf5 format supported by FastReg()
#' @param dataFile Path of the input text file.
#' @param h5Dir Path of a new directory to write generated hdf5 files.
#' @param headerRow Default: 1. Row number of header line.
#' @param idCol Default 1. Column number that contains identifiers.
#' @param dataCol Default: 2. Column number that contains first data entry.
#' @param buffSize Default: 1.0. Size in Gb of memory available for data conversion.
#' @param transpose Default: FALSE. Boolean value that indicates whether text data should be transposed. Set TRUE when columns represent individuals and rows represent POIs.
#' @param chunkEdge Default: 100. Specify size of data chunks to be applied to output hdf5 files. Must be less than or equal to values of poiPerFile and greater than zero.
#' @param vcf Default: FALSE. Indicate whether input file is in vcf format. If TRUE, headerRow, idCol, dataCol and transpose values are ignored.
#' @param delimiter Specify a string of all delimiters used in the text file to separate columns. Default will split entries by space or tab.
#' @param gz Default: FALSE. Indicate whether input text file is gzip-compressed.
#' @param poiPerFile Default: -1. Indicate the number of POIs to write to each output hdf5 file. A value of -1 indicates the count should be calculated based on available system resources.
#' @param singleFile Default: FALSE. Indicate whether a single hdf5 file should be produced rather than a series.
#' @param serverThreads Default: -1. Indicate the total number of CPU threads available for FastReg - utilized to determine the optimal number of output files. A value of -1 indicates the count should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.
#' @param serverMem Default: -1.0. Indicate the total memory in Gb available for FastReg - utilized to determine the optimal number of output files. A value of -1.0 indicates the total should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.
#' @return 0 = success, 1 = failure
#' @export

FastRegImport <- function(
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
  FastRegImportCpp(
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
