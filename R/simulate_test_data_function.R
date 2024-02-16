#' simulate_test_dataset function to generate test dataset for evaluation of FastReg
#' @param num.poi an integer denoting number of predictors of interest (SNP calls)
#' @param num.ind an integer denoting number of individuals
#' @param seed an integer denoting random number generate seed
#' @param coeff.sd numeric value denoting the standard deviation parameter for regression coefficients for covariates (default=0, regression coefficients for covarites are set to zero)
#' @param bin.resp.mean numeric value between 0 and 1 denoting baseline incidence rate for binary response
#' @param num.resp.mean numeric value denoting overall mean of numeric response
#' @param num.resp.sd numeric value denoting overall standard deviation of numeric response
#' @param poi.type character must be either "genotype" or "dosage"
#' @param poi.chunk.size an integer denoting poi chunk size used during H5 file generation
#' @param poi.compression.level an integer denoting compression level used during H5 file creation
#' @param data.dir output directory
#' @param prefix prefix for all files created
#' @param poi.file.type 'h5' or 'txt'
#' @param verbose logical to control display of progress messages (default=TRUE)
#' @import stats
#' @import parallel
#' @import data.table
#' @import rhdf5
#' @return a list consisting of dataset size (num.poi, num.ind) as well as regression coefficients used for both binary and numeric response
simulate_test_dataset <- function(num.poi = 50000,
                                  num.ind = 5000,
                                  seed = 12133,
                                  coeff.sd = 0,
                                  bin.resp.mean = 0.2,
                                  num.resp.mean = 24,
                                  num.resp.sd = 5,
                                  poi.type = "genotype",
                                  poi.chunk.size = 100,
                                  poi.compression.level = 7,
                                  data.dir = ".",
                                  data.num.chunks = 10,
                                  prefix = "testdata_5k_by_50k",
                                  verbose = TRUE,
                                  poi.file.type = "txt") {
  # num.poi <- 50000
  # num.ind <- 5000
  covariates <- c("age", "sex", "eth", "treatment", "severity")

  # seed <- 12133;
  # poi.type <- "genotypes"
  # poi.chunk.size <- 100
  # poi.compression.level <- 7
  # data.dir <- "../input/"
  # prefix = "testdata_5k_by_50k"

  dir.create(data.dir, showWarnings = FALSE)

  covar.file <- file.path(data.dir, paste0(prefix, ".covar.txt"))
  num.pheno.file <- file.path(data.dir, paste0(prefix, ".num.pheno.txt"))
  bin.pheno.file <- file.path(data.dir, paste0(prefix, ".bin.pheno.txt"))
  poi.txt.file <- file.path(data.dir, paste0(prefix, ".poi.txt"))
  poi.data.dir <- paste0(data.dir, "/h5")
  dir.create(poi.data.dir, showWarnings = FALSE)
  poi.file <- file.path(poi.data.dir, paste0(c("/", prefix, ".poi.h5")))
  poi.subset.file <- file.path(data.dir, paste0(prefix, ".poi.subset.txt"))
  subject.subset.file <- file.path(data.dir, paste0(prefix, ".sample.subset.txt"))


  set.seed(seed)

  ind.id <- paste0("IND", formatC(1:num.ind, format = "d", flag = "0", digits = floor(log10(num.ind))))
  poi.id <- paste0("rs", formatC(1:num.poi, format = "d", flag = "0", digits = floor(log10(num.poi))))


  write.table(data.frame("ind" = sort(sample(ind.id, size = 5, replace = FALSE)), stringsAsFactors = FALSE),
    file = subject.subset.file, quote = FALSE, na = "", row.names = FALSE, col.names = TRUE
  )
  write.table(data.frame("rsID" = sort(sample(poi.id, size = 50, replace = FALSE)), stringsAsFactors = FALSE),
    file = poi.subset.file, quote = FALSE, na = "", row.names = FALSE, col.names = TRUE
  )


  covar.df <- data.frame(
    "ind" = ind.id, "age" = as.integer(rnorm(num.ind, 50, 7)),
    "sex" = c("M", "F")[ceiling(2 * runif(num.ind))],
    "eth" = c("eur", "afr", "asi")[ceiling(3 * runif(num.ind))],
    "treatment" = c("Test", "Placebo")[ceiling(2 * runif(num.ind))],
    "severity" = c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme")[ceiling(6 * runif(num.ind))],
    stringsAsFactors = FALSE, check.names = FALSE
  )

  write.table(covar.df, file = covar.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  remove(list = "covar.df")
  if (verbose) cat("generated covariate file\n")

  num.resp.beta0 <- num.resp.mean
  bin.resp.beta0 <- log(bin.resp.mean / (1 - bin.resp.mean))


  if (coeff.sd > 0) {
    X <- cbind(
      "Age" = (covar.df$Age - mean(covar.df$Age)) / sd(covar.df$Age),
      "sex:Male vs. Female" = 1 * (covar.df$sex == "M"),
      "eth:afr vs. eur" = 1 * (covar.df$eth == "afr"), "eth:asi vs. eur" = 1 * (covar.df$eth == "asi"),
      "treatment:Test vs. Placebo" = 1 * (covar.df$treatment == "Test"),
      "severity:Low vs. Very Low" = 1 * (covar.df$severity == "Low"),
      "severity:Moderate vs. Very Low" = 1 * (covar.df$severity == "Moderate"),
      "severity:High vs. Very Low" = 1 * (covar.df$severity == "High"),
      "severity:Very High vs. Very Low" = 1 * (covar.df$severity == "Very High"),
      "severity:Extreme vs. Very Low" = 1 * (covar.df$severity == "Extreme")
    )


    num_resp.beta <- matrix(rnorm(ncol(X), mean = 0, sd = coeff.sd), ncol = 1, dimnames = list(colnames(X), NULL))
    bin_resp.beta <- matrix(rnorm(ncol(X), mean = 0, sd = coeff.sd), ncol = 1, dimnames = list(colnames(X), NULL))
    num.resp <- rnorm(num.ind, mean = num.resp.beta0 + X %*% num_resp.beta, sd = num.resp.sd)
    bin.resp <- rbinom(num.ind, prob = 1 / (1 + exp(-(bin.resp.beta0 + X %*% bin_resp.beta))), size = 1)
    remove(list = "X")
  } else {
    num.resp <- rnorm(num.ind, mean = num.resp.beta0, sd = num.resp.sd)
    bin.resp <- rbinom(num.ind, prob = 1 / (1 + exp(-(bin.resp.beta0))), size = 1)
    num_resp.beta <- NULL
    bin_resp.beta <- NULL
  }

  pheno.df <- data.frame("ind" = ind.id, "num.resp" = num.resp, "bin.resp" = bin.resp)


  write.table(pheno.df[, c(1, 2)], file = num.pheno.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  write.table(pheno.df[, c(1, 3)], file = bin.pheno.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  remove(list = "pheno.df")
  if (verbose) cat("generated phenotype files\n")
  num.poi.blocks <- ceiling(num.poi / poi.chunk.size)
  if (poi.file.type == "h5") {
    write.h5(
      data.num.chunks,
      prefix,
      verbose,
      ind.id,
      poi.id,
      poi.data.dir,
      poi.compression.level,
      poi.chunk.size,
      poi.type,
      num.ind,
      num.poi
    )
  } else {
    values <- generate_values(num.ind, 1:num.poi, poi.type)
    values <- t(values)
    colnames(values) <- c("", ind.id)
    write.table(values, file = poi.txt.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  }

  # file.con[["values"]] <- t(values);
  # file.con$close_all();

  invisible(list("num.poi" = num.poi, "num.ind" = num.ind, "num.resp.beta0" = num.resp.beta0, "num_resp.beta" = num_resp.beta, "num.resp.sd" = num.resp.sd, "bin.resp.beta0" = bin.resp.beta0, "bin_resp.beta" = bin_resp.beta))
}

generate_values <- function(num.ind, poi.indices, poi.type) {
  block.index <- poi.indices
  block.size <- length(block.index)
  # miss.rate <- runif(block.size, min = 0, max = 0.1)

  if (poi.type == "genotype") {
    values <- matrix(0, ncol = block.size, nrow = num.ind)
  } else {
    values <- matrix(0.000, ncol = block.size, nrow = num.ind)
  }
  maf <- runif(block.size, min = 0.05, max = 0.5)
  miss.rate <- runif(block.size, min = 0, max = 0.1)
  for (i in 1:block.size) {
    if (poi.type == "genotype") {
      dosage.val <- runif(num.ind)
      geno.val <- integer(num.ind)
      geno.val[dosage.val < (1 - maf[i])^2] <- 0
      geno.val[((1 - maf[i])^2 < dosage.val) & (dosage.val < (1 - maf[i]^2))] <- 1
      geno.val[dosage.val > (1 - maf[i]^2)] <- 2
      geno.val[runif(num.ind)>1-miss.rate[i]] <- NA;
      values[, i] <- geno.val
    } else {
      geno.val <- ((1 - maf[i])^2) * runif(num.ind, min = 0, max = 0.4) + 2 * (1 - maf[i]) * maf[i] * runif(num.ind, min = 0.25, max = 0.75) + (maf[i]^2) * runif(num.ind, min = 0.6, max = 1)
      geno.val[runif(num.ind) > 1 - miss.rate[i]] <- NA
      values[, i] <- geno.val
      # miss.rate <- runif(block.size, min = 0, max = 0.1)
    }
  }
  values
}

write.h5 <- function(
    data.num.chunks,
    prefix,
    verbose,
    ind.id,
    poi.id,
    poi.data.dir,
    poi.compression.level,
    poi.chunk.size,
    poi.type,
    num.ind,
    num.poi) {
  total.pois.per.file <- ceiling(num.poi / data.num.chunks)
  pois.per.file <- rep(total.pois.per.file, data.num.chunks)
  pois.per.file[data.num.chunks] <- num.poi - total.pois.per.file * (data.num.chunks - 1)
  num.poi.chunk <- ceiling(num.poi / data.num.chunks)
  for (k in 1:data.num.chunks) {
    if (k == data.num.chunks) {
        # For the last chunk, adjust the end.index to ensure it doesn't exceed the total number of POIs
        start.index <- ((k - 1) * total.pois.per.file) + 1
        end.index <- num.poi  # Ensure the last chunk includes all remaining POIs
    } else {
        start.index <- ((k - 1) * total.pois.per.file) + 1
        end.index <- k * total.pois.per.file
    }
    # start.index <- sum(pois.per.file[1:(k - 1)]) + 1
    # end.index <- sum(pois.per.file[1:k])

    poi.file <- file.path(poi.data.dir, paste0(prefix, ".", k, ".poi.h5"))
    if (file.exists(poi.file)) unlink(poi.file)
    # file.con <- H5File$new(poi.file, mode="w");
    h5createFile(file = poi.file)

    if (verbose) cat("added ind id into h5 file\n")
    # file.con[["individuals"]] <- ind.id
    # h5createGroup(file=poi.file, "individuals");
    h5write(ind.id, level = 4, file = poi.file, name = "individuals")

    if (verbose) cat("added poi id into h5 file\n")

    # file.con[["predictors_of_interest"]] <- poi.id
    #  h5createGroup(file=poi.file, "predictors_of_interest");
    h5write(poi.id[start.index:end.index], level = 4, file = poi.file, name = "predictors_of_interest")

    h5createDataset(
      file = poi.file,
      dataset = "values",
      dims = c(num.ind, num.poi.chunk),
      storage.mode = ifelse(poi.type == "genotype", "integer", "double"),
      chunk = c(num.ind, min(poi.chunk.size, pois.per.file[k])),
      level = poi.compression.level
    )
    for (block.start in seq(start.index, end.index, by = poi.chunk.size)) {
      block.end <- min(block.start + poi.chunk.size - 1, end.index)
      values <- generate_values(num.ind, block.start:block.end, poi.type)
      h5write(values, file = poi.file, name = "values", start = c(1, block.start - start.index + 1), count = c(num.ind, block.end - block.start + 1))
    }
  }
  if (verbose) cat("added POI values into h5 file\n")
  h5closeAll()
}
