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
#' @param max.openmp.threads Default: 1. Overrides number of threads used. By default FastReg will use 1 openmp threads.
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
#' @param max.blas.threads Default: -1. Set the number of BLAS threads that Eigen/Armadillo will use. The default value will use the system/BLAS library default.
#' @param max.workers Default: 0. Number of processes FastReg will spawn. By default FastReg will spawn min(# of poi files, # of available workers)
#' @return numeric denoting elapsed.time
#' @export
#' @import rhdf5
#' @import parallel
VLA <- function(
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
    max.openmp.threads = 1,
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
    max.blas.threads = -1,
    max.workers = 0) {
  if (max.openmp.threads <= 0) {
    cat("Error: max.openmp.threads must be a positive integer.\n")
    return(FALSE)
  }

  # if (max.blas.threads > get_num_procs()) {
  #   cat("Error: max.blas.threads cannot be greater than number of avaialable processors.\n")
  #   return(FALSE)
  # }

  # if (max.blas.threads > 0) {
  #   blas_set_num_threads(max.blas.threads)
  # }

  if (max.workers < 0) {
    cat("Error: max.workers must be a positive integer.\n")
    return(FALSE)
  }

  # FastRegVLA(
  #   phenotype,
  #   regression.type,
  #   Pvalue.dist,
  #   output.exclude.covar,
  #   maf.threshold,
  #   hwe.threshold,
  #   no.intercept,
  #   colinearity.rsq,
  #   poi.block.size,
  #   max.iter,
  #   rel.conv.tolerance,
  #   abs.conv.tolerance,
  #   max.openmp.threads,
  #   pheno.file,
  #   pheno.rowname.cols,
  #   pheno.file.delim,
  #   covar.file,
  #   covar.rowname.cols,
  #   covar.file.delim,
  #   POI.file.dir,
  #   POI.file.delim,
  #   POI.file.format,
  #   POI.type,
  #   POI.effect.type,
  #   covariates,
  #   covariate.type,
  #   covariate.standardize,
  #   covariate.levels,
  #   covariate.ref.level,
  #   POI.covar.interactions,
  #   split.by,
  #   output.dir,
  #   compress.results,
  #   max.workers
  # )
}

#' TextToH5 a function to convert textual data to hdf5 format supported by FastReg()
#' @param data.file Path of the input text file.
#' @param h5.dir Path of a new directory to write generated hdf5 files.
#' @param header.row Default: 1. Row number of header line.
#' @param id.col Default 1. Column number that contains identifiers.
#' @param data.col Default: 2. Column number that contains first data entry.
#' @param buff.size Default: 1.0. Size in Gb of memory available for data conversion.
#' @param transpose Default: FALSE. Boolean value that indicates whether text data should be transposed. Set TRUE when columns represent individuals and rows represent POIs.
#' @param chunk.edge Default: 100. Specify size of data chunks to be applied to output hdf5 files. Must be less than or equal to values of poi.per.file and greater than zero.
#' @param vcf Default: FALSE. Indicate whether input file is in vcf format. If TRUE, header.row, id.col, data.col and transpose values are ignored.
#' @param delimiter Specify a string of all delimiters used in the text file to separate columns. Default will split entries by space or tab.
#' @param gz Default: FALSE. Indicate whether input text file is gzip-compressed.
#' @param poi.per.file Default: -1. Indicate the number of POIs to write to each output hdf5 file. A value of -1 indicates the count should be calculated based on available system resources.
#' @param single.file Default: FALSE. Indicate whether a single hdf5 file should be produced rather than a series.
#' @param server.threads Default: -1. Indicate the total number of CPU threads available for FastReg - utilized to determine the optimal number of output files. A value of -1 indicates the count should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.
#' @param server.mem Default: -1.0. Indicate the total memory in Gb available for FastReg - utilized to determine the optimal number of output files. A value of -1.0 indicates the total should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.
#' @return 0 = success, 1 = failure
#' @export

TextToH5 <- function(
    data.file,
    h5.dir,
    header.row = 1,
    id.col = 1,
    data.col = 2,
    buff.size = 1.0,
    transpose = FALSE,
    chunk.edge = 100,
    vcf = FALSE,
    delimiter = " \t",
    gz = FALSE,
    poi.per.file = -1,
    single.file = FALSE,
    server.threads = -1,
    server.mem = -1.0) {
  if (!is.character(data.file) || length(data.file) > 1) {
    cat("Error: data.file must be a single-member character vector\n")
    return(FALSE)
  }
  if (!is.character(h5.dir) || length(h5.dir) > 1) {
    cat("Error: h5.dir must be a single-member character vector\n")
    return(FALSE)
  } else if (file.exists(h5.dir)) {
    cat("Error: file or folder named \"", h5.dir, "\" already exists\n", sep = "")
    return(FALSE)
  }
  dir.create(h5.dir)
  if (!is.numeric(header.row) || header.row != as.integer(header.row) || length(header.row) > 1 || header.row < 1) {
    cat("Error: header.row must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  header.row <- as.integer(header.row)
  if (!is.numeric(id.col) || id.col != as.integer(id.col) || length(id.col) > 1 || id.col < 1) {
    cat("Error: id.col must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  id.col <- as.integer(id.col)
  if (!is.numeric(data.col) || data.col != as.integer(data.col) || length(data.col) > 1 || data.col <= id.col) {
    cat("Error: header.row must be a single-member integer vector with value > id.col\n")
    return(FALSE)
  }
  data.col <- as.integer(data.col)
  if (!is.numeric(buff.size) || length(buff.size) > 1 || buff.size <= 0) {
    cat("Error: buff.size must be a single-member numeric vector with value > 0\n")
    return(FALSE)
  }
  buff.size <- as.double(buff.size)
  if (!is.logical(transpose) || length(transpose) > 1) {
    cat("Error: transpose must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(chunk.edge) || chunk.edge != as.integer(chunk.edge) || length(chunk.edge) > 1 || chunk.edge < 1) {
    cat("Error: chunk.edge must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  chunk.edge <- as.integer(chunk.edge)
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
  if (!is.numeric(poi.per.file) || poi.per.file != as.integer(poi.per.file) || length(poi.per.file) > 1 || poi.per.file < chunk.edge) {
    if (poi.per.file != -1) {
      cat("Error: poi.per.file must be a single-member integer vector with value >= chunk.edge\n")
      return(FALSE)
    }
  }
  poi.per.file <- as.integer(poi.per.file)
  if (!is.logical(single.file) || length(single.file) > 1) {
    cat("Error: single.file must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(server.threads) || server.threads != as.integer(server.threads) || length(server.threads) > 1) {
    cat("Error: server.threads must be a single-member integer vector\n")
    return(FALSE)
  }
  server.threads <- as.integer(server.threads)
  if (!is.numeric(server.mem) || length(server.mem) > 1) {
    cat("Error: server.mem must be a single-member numeric vector\n")
    return(FALSE)
  }
  server.mem <- as.double(server.mem)
  # FastRegImportCpp(
  #   data.file,
  #   h5.dir,
  #   header.row,
  #   id.col,
  #   data.col,
  #   buff.size,
  #   transpose,
  #   chunk.edge,
  #   vcf,
  #   delimiter,
  #   gz,
  #   poi.per.file,
  #   single.file,
  #   server.threads,
  #   server.mem
  # )
}
