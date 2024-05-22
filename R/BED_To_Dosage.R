#' Bed_To_dosage function to convert previously generated H5 files into plink files
#' @param input.dir folder location where input files reside
#' @param output.dir folder location where plink files will be saved
#' @param prefix a character denoting which set of h5 files to convert (first set is used by default)
#' @param poi.chunk.size integer denoting chunk.size to be used during conversion
#' @return integer 1 when successful
#' @import BEDMatrix
#' @import data.table

Bed_To_dosage <- function(input.dir = ".", output.dir = NULL, prefix = NULL, poi.chunk.size = 1000) {
  
  if (is.null(output.dir)) output.dir <- input.dir
  
  available.files <- list.files(path = input.dir, pattern = "[.]bed")
  if (length(available.files) == 0) stop("input.dir does not consist of any .bed file")
  
  if (is.null(prefix)) prefix <- sub("[.]bed", "", basename(available.files)[1])
  if (length(prefix) > 1) {
    warning("only first prefix used")
    prefix <- prefix[1]
  }
  
  src.files <- file.path(input.dir, paste0(prefix, ".", c("bim", "fam", "bed")))
  out.files <- file.path(output.dir, paste0(prefix, ".", "dosage.txt"))

  fam_data <- fread(src.files[2], header = FALSE)
  setnames(fam_data, c("FID", "IID", "PID", "MID", "Sex", "Phenotype"))
  fid_headers <- fam_data$FID
  iid_header <- fam_data$IID
  print(head(fam_data))
  num.poi <- as.integer(strsplit(system(command = paste("wc -l", src.files[1]), intern = TRUE), split = " ")[[1]][1])
  num.ind <- nrow(fam_data)
  
  num.poi.chunks <- ceiling(num.poi / poi.chunk.size)
  poi.chunk.size <- ceiling(num.poi / num.poi.chunks)
  
  bm <- BEDMatrix(src.files[3], n = num.ind, p = num.poi)
  sample_headers <- unlist(mapply(function(f, i) c(f, i), fam_data$FID, fam_data$IID, SIMPLIFY = FALSE))
  for (i in 1:num.poi.chunks) {
    start.pos <- (i - 1) * poi.chunk.size + 1
    end.pos <- min(i * poi.chunk.size, num.poi)
    nr <- end.pos - start.pos + 1
    snp_info <- fread(src.files[1], skip = start.pos - 1, nrows = nr, header = FALSE)
    setnames(snp_info, c("Chromosome", "SNP", "GenDist", "Position", "A1", "A2"))
    
    poi_block <- t(bm[, seq(start.pos, end.pos), drop = FALSE])

    dosage_block <- matrix("NA", nrow = nr, ncol = num.ind * 2)
    for (genotype in c("0", "1", "2")) {
        genotype_indices <- poi_block == genotype
        dosage_block[genotype_indices] <- ifelse(genotype == "0", "0.00", ifelse(genotype == "1", "1.00", "2.00"))
    }
    dosage_data <- cbind(snp_info[, .(SNP, A1, A2)], as.data.table(dosage_block))
    expected_colnames <- c("SNP", "A1", "A2", sample_headers) 
    if (ncol(dosage_data) != length(expected_colnames)) {
      stop("Number of columns in the data table does not match the number of names assigned.")
    }
    setnames(dosage_data, expected_colnames)

    # Write to file
    fwrite(dosage_data, file = out.files, append = (i != 1), sep = "\t", col.names = (i == 1))
  }
  
  invisible(1)
}
