#' perform.GLM.test function to validate FastReg results
#' @param num.poi an integer denoting number of POI to be evaluated
#' @param prefix a character string denoting filename prefix (same as one used by simulate_test_dataset function)
#' @param input.dir a character string denoting folder location input files
#' @param output.dir a character string denoting folder location of output files;
#' @param transpose.h5 a logical denoting whether to transpose POI matrix
#' @param na.values vector of POI values that should treated as missing
#' @param mc.cores an integer denoting number of cores utilzed during parallel processing
#' @return invisible TRUE when successful (POI level output for logistic and linear regression are saved in separate tsv files)
#' @export
#' @import parallel
perform.GLM.test <- function(num.poi = 10, prefix = "testdata_500_by_500", input.dir = ".", output.dir = ".", transpose.h5 = FALSE, na.values = NA, mc.cores = 1) {
  x.df <- read.delim(file = file.path(input.dir, paste0(prefix, ".covar.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  yb.df <- read.delim(file = file.path(input.dir, paste0(prefix, ".bin.pheno.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  yn.df <- read.delim(file = file.path(input.dir, paste0(prefix, ".num.pheno.txt")), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  if (!transpose.h5) {
    z.df <- h5read(file = file.path(input.dir, paste0(prefix, ".poi.h5")), "/values", index = list(NULL, 1:num.poi))
  } else {
    z.df <- t(h5read(file = file.path(input.dir, paste0(prefix, ".poi.h5")), "/values", index = list(1:num.poi, NULL)))
  }

  z.df[z.df %in% na.values] <- NA


  dimnames(z.df) <- list(
    h5read(file = file.path(input.dir, paste0(prefix, ".poi.h5")), "/individuals"),
    h5read(file = file.path(input.dir, paste0(prefix, ".poi.h5")), "/predictors_of_interest")[1:num.poi]
  )

  z.df <- data.frame("ind" = rownames(z.df), z.df, stringsAsFactors = FALSE, check.names = FALSE)

  df <- merge(yb.df, yn.df, by = "ind")
  df <- merge(df, x.df, by = "ind")
  df <- merge(df, z.df, by = "ind")

  df$age <- (df$age - mean(df$age, na.rm = TRUE)) / sd(df$age, na.rm = TRUE)
  df$sex <- factor(df$sex, levels = c("F", "M"), ordered = FALSE)
  df$eth <- factor(df$eth, levels = c("eur", "afr", "asi"), ordered = FALSE)
  df$treatment <- factor(df$treatment, levels = c("Placebo", "Test"), ordered = FALSE)
  df$severity <- factor(df$severity, levels = c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme"), ordered = FALSE)


  poi.names <- setdiff(colnames(z.df), "ind")

  remove(list = c("x.df", "yn.df", "yb.df", "z.df"))

  # cat(poi.names, "\n")
  process.poi <- function(poi.name, type) {
    dat <- df[, c(type, "age", "sex", "eth", "treatment", "severity", poi.name)]
    colnames(dat)[ncol(dat)] <- "POI"

    family <- ifelse(type == "bin.resp", "binomial", "gaussian")
    pval <- ifelse(type == "bin.resp", "Pr(>|z|)", "Pr(>|t|)")
    keep.ind <- which(rowSums(do.call(cbind, lapply(dat[c(type, "age", "sex", "eth", "treatment", "severity", "POI")], FUN = is.na))) == 0)


    glm.obj <- glm(as.formula(paste0(type, "~age+sex+eth+treatment+severity+POI")), data = dat[keep.ind, ], family = family)

    glm.summary.obj <- summary(glm.obj)
    k <- nrow(glm.summary.obj$coef)

    res <- data.frame(
      poi = rep(poi.name, k),
      DF = rep(glm.summary.obj$df.residual, k),
      par = rownames(glm.summary.obj$coef),
      glm.summary.obj$coeff[, c("Estimate", "Std. Error", pval)],
      check.names = FALSE, stringsAsFactors = FALSE
    )
    res$"NegLog10 P-value" <- -(pt(abs(res$Estimate / res$"Std. Error"), df = res$DF, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
    res$"NegLog10 P-value (indirect)" <- -log(2 * (1 - pt(abs(res$Estimate / res$"Std. Error"), df = res$DF))) / log(10)


    return(res)
  }

  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)


  for (type in c("bin.resp", "num.resp")) {
    out <- do.call(rbind, mclapply(poi.names, FUN = process.poi, mc.cores = mc.cores, mc.preschedule = TRUE, type = type))
    write.table(out, file = file.path(output.dir, paste0(prefix, ".", type, "-glm-Output.txt")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
  }


  invisible(TRUE)
}
