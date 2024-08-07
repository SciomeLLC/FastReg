#' H5_To_Bed function to convert previously generated H5 files into plink files
#' @param input.dir folder location where input files reside
#' @param output.dir folder location where plink files will be saved
#' @param prefix a character denoting which set of h5 files to convert (first set is used by default)
#' @param poi.chunk.size integer denoting chunk.size to be used during conversion
#' @param seed random number generator seed
#' @param num.poi number of poi used during conversion (default=NULL, equivalent to all)
#' @param num.ind number of ind used during conversion (default=NULL, equivalent to all)
#' @return integer 1 when successful
#' @import stats
#' @import parallel
H5_To_Bed <- function(input.dir = ".", output.dir = NULL, prefix = NULL, poi.chunk.size = 1000, seed = 12133, num.poi = NULL, num.ind = NULL) {
  if (is.null(output.dir)) output.dir <- input.dir

  available.files <- list.files(path = input.dir, pattern = "[.]poi[.]h5")

  if (length(available.files) == 0) stop("input.dir does not consist of any .poi.h5 file")
  if (is.null(prefix)) prefix <- sub("[.]poi[.]h5", "", basename(available.files)[1])

  if (length(prefix) > 1) {
    warning("only first prefix used")
    prefix <- prefix[1]
  }

  src.files <- file.path(input.dir, paste0(prefix, ".", c("poi.h5", "covar.plink.txt", "num.pheno.plink.txt", "bin.pheno.plink.txt")))
  out.files <- file.path(output.dir, paste0(prefix, ".", c("bim", "fam", "bed")))

  if (any(!file.exists(src.files))) stop("critical files are missing")

  ### Hg38 chromsome lengths
  seqs <- c(
    "chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559, "chr4" = 190214555, "chr5" = 181538259,
    "chr6" = 170805979, "chr7" = 159345973, "chr8" = 145138636, "chr9" = 138394717, "chr10" = 133797422,
    "chr11" = 135086622, "chr12" = 133275309, "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189,
    "chr16" = 90338345, "chr17" = 83257441, "chr18" = 80373285, "chr19" = 58617616, "chr20" = 64444167,
    "chr21" = 46709983, "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415, "chrMT" = 16569
  )

  seq.names <- names(seqs)
  seq.names <- gsub("chr", "", seq.names)
  num.seq <- length(seqs)

  poi.id <- h5read(file = src.files[1], "predictors_of_interest")
  ind.id <- h5read(file = src.files[1], "individuals")

  if (!is.null(num.ind)) ind.id <- ind.id[1:num.ind]
  if (!is.null(num.poi)) poi.id <- poi.id[1:num.poi]

  num.poi <- length(poi.id)
  num.ind <- length(ind.id)

  set.seed(seed)

  cat("Generating .bim file\n")
  poi.seq <- as.integer(cut(sort(runif(num.poi)), breaks = c(0, cumsum(as.numeric(seqs / sum(as.numeric(seqs)))))))
  poi.pos <- integer(num.poi)
  poi.ref <- poi.alt <- character(num.poi)
  for (seq in 1:num.seq) {
    w <- which(poi.seq == seq)
    num.w <- length(w)
    if (num.w > 0) {
      poi.pos[w] <- sort((unique(ceiling(seqs[seq] * runif(num.w * 2))))[1:num.w])
      rf <- ceiling(4 * runif(num.w))
      poi.ref[w] <- c("A", "T", "C", "G")[rf]
      poi.alt[w] <- c("A", "T", "C", "G")[4 - rf + 1]
    }
  }

  poi.df <- as.data.table(data.frame(CHROM = seq.names[poi.seq], ID = poi.id, POS = integer(num.poi), "COORD" = poi.pos, "ALT" = poi.alt, "REF" = poi.ref, stringsAsFactors = FALSE, check.names = FALSE))
  remove(list = c("poi.seq", "poi.pos", "poi.ref", "poi.alt", "w", "num.w", "rf", "seq", "seqs", "num.seq", "seq.names"))
  write.table(poi.df, file = out.files[1], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "")
  remove(list = "poi.df")

  cat("Generating .fam file\n")
  ## load covar file
  covar <- as.data.frame(fread(file = src.files[2], sep = "\t", header = TRUE), as.is = TRUE)
  pheno.num <- as.data.frame(fread(file = src.files[3], sep = "\t", header = TRUE), as.is = TRUE)
  pheno.bin <- as.data.frame(fread(file = src.files[4], sep = "\t", header = TRUE), as.is = TRUE)
  
  # pheno.bin <- rename(pheno.bin, IID = ind)
  # pheno.num <- rename(pheno.num, IID = ind)
  covar <- merge(merge(covar, pheno.num, by = "IID", all = TRUE), pheno.bin, by = "IID", all = TRUE)

  fam <- data.frame(
    `#FID` = covar$`#FID`, 
    IID = ind.id, 
    WID_F = integer(num.ind), 
    WID_M = integer(num.ind), 
    SEX = integer(num.ind), 
    PHENO = integer(num.ind), 
    check.names = FALSE
  )

  # common.ind <- which(fam$IID %in% as.character(covar$IID))
  common.ind <- which(fam$IID %in% covar$IID)
  if (length(common.ind) == 0 || is.null(common.ind)) stop("No matching IDs found between fam and covar dataframes.")

  # ind.index <- match(fam$`#FID`[common.ind], as.character(covar$`#FID`))
  ind.index <- match(fam$`#FID`[common.ind], covar$`#FID`)
  sex_assignment <- (1 * as.integer(covar[ind.index, "sex"] == "M")) + (2 * as.integer(covar[ind.index, "sex"] == "F"))

  fam$SEX[common.ind] <- sex_assignment
  fam$SEX[is.na(fam$SEX)] <- 0
  fam$PHENO[is.na(fam$SEX)] <- 0
  # Check if any replacements have length zero
  if (length(fam$SEX[common.ind]) == 0 || length(fam$PHENO[common.ind]) == 0) {
    stop("Replacement has length zero due to incorrect indexing or data misalignment.")
  }
  write.table(fam, file = out.files[2], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "")
  remove(list = "fam")

  cat("Generating .bed file\n")
  if (file.exists(out.files[3])) unlink(out.files[3])
  con <- file(out.files[3], "wb")
  poi.chunks <- ceiling(num.poi / poi.chunk.size)

  nBytes <- ceiling(num.ind / 4)

  d <- matrix(c(
    0, 1, 1, 0, 1, 1, 0, 0,
    0, 0, 0, 1, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1
  ), nrow = 8)

  xd <- as.raw(drop(matrix(2^(7:0), nrow = 1, ncol = 8) %*% d))
  temp_output <- writeBin(xd, con = con)

  base.vec <- matrix(2^(0:7), nrow = 1, ncol = 8)

  for (j in 1:poi.chunks) {
    poi.index <- ((j - 1) * poi.chunk.size + 1):(min(j * poi.chunk.size, num.poi))
    chunk.size <- length(poi.index)
    poi.chunk <- h5read(file = src.files[1], "values", index = list(1:num.ind, poi.index))

    xd <- matrix(0, nrow = nBytes, ncol = chunk.size)

    for (k in 1:chunk.size) {
      d <- matrix(0, nrow = 8, ncol = nBytes)
      ip <- which(poi.chunk[, k] == 1)
      if (length(ip) > 0) {
        d[2 * (ip - 1) + 1] <- 0
        d[2 * (ip - 1) + 2] <- 1
      }
      ip <- which(poi.chunk[, k] == 0)
      if (length(ip) > 0) {
        d[2 * (ip - 1) + 1] <- 1
        d[2 * (ip - 1) + 2] <- 1
      }
      ip <- which(is.na(poi.chunk[, k]))
      if (length(ip) > 0) {
        d[2 * (ip - 1) + 1] <- 1
        d[2 * (ip - 1) + 2] <- 0
      }
      xd[, k] <- drop(base.vec %*% d)
      remove(list = "d")
    }

    # print(xd[, k])
    data <- as.raw(c(xd))
    # print(data)
    writeBin(data, con=con); flush(con);

    remove(list = c("poi.chunk", "xd"))
  }

  close(con)

  return(invisible(1))
}