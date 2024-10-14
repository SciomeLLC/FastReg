#' create_test_dataset function to generate test dataset for evaluation of FastReg
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
#' @param data.num.chunks number of files to split the data into
#' @param prefix prefix for all files created
#' @param covariates list of covariates to generate. Supported types are: age, sex, treatment, eth, severity
#' @param poi.file.type 'h5' or 'txt'
#' @param verbose logical to control display of progress messages (default=TRUE)
#' @return a list consisting of dataset size (num.poi, num.ind) as well as regression coefficients used for both binary and numeric response
#' @import stats
#' @import parallel
#' @import rhdf5
#' @import utils
#' @import doParallel
#' @import foreach
create_test_dataset <- function(num.poi = 50000,
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
                                  covariates = c("age", "sex", "treatment"),
                                  poi.file.type = "txt",
                                  verbose = TRUE) {
  dir.create(data.dir, showWarnings = FALSE)

  covar.file <- file.path(data.dir, paste0(prefix, ".covar.txt"))
  covar.plink.file <- file.path(data.dir, paste0(prefix, ".covar.plink.txt"))
  num.pheno.file <- file.path(data.dir, paste0(prefix, ".num.pheno.txt"))
  bin.pheno.plink.file <- file.path(data.dir, paste0(prefix, ".bin.pheno.plink.txt"))
  num.pheno.plink.file <- file.path(data.dir, paste0(prefix, ".num.pheno.plink.txt"))
  bin.pheno.file <- file.path(data.dir, paste0(prefix, ".bin.pheno.txt"))
  poi.txt.file <- file.path(data.dir, paste0(prefix, ".poi.txt"))
  poi.data.dir <- paste0(data.dir, "/h5")
  poi.file <- file.path(poi.data.dir, paste0(c("/", prefix, ".poi.h5")))
  poi.subset.file <- file.path(data.dir, paste0(prefix, ".poi.subset.txt"))
  subject.subset.file <- file.path(data.dir, paste0(prefix, ".sample.subset.txt"))


  set.seed(seed)
  fid.id <- paste0("FID", formatC(1:num.ind, format = "d", flag = "0", digits = floor(log10(num.ind))))
  ind.id <- paste0("IND", formatC(1:num.ind, format = "d", flag = "0", digits = floor(log10(num.ind))))
  poi.id <- paste0("rs", formatC(1:num.poi, format = "d", flag = "0", digits = floor(log10(num.poi))))


  write.table(data.frame("IID" = sort(sample(ind.id, size = 5, replace = FALSE)), stringsAsFactors = FALSE),
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
  covar.df <- covar.df[, c("ind", covariates)]
  write.table(covar.df, file = covar.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  if (verbose) cat("generated covariate file\n")

  names(covar.df)[names(covar.df) == "ind"] <- "IID"
  covar.plink.df <- data.frame("#FID" = fid.id, covar.df, check.names = FALSE)
  # print(head(covar.plink.df))
  write.table(covar.plink.df, file = covar.plink.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  if (verbose) cat("generated covariate file for plink\n")
  remove(list = "covar.df")
  remove(list = "covar.plink.df")

  num.resp.beta0 <- num.resp.mean
  bin.resp.beta0 <- log(bin.resp.mean / (1 - bin.resp.mean))

  if (coeff.sd > 0) {
    X <- matrix(nrow = num.ind, ncol = 0);
    if ("age" %in% covariates) {
      X <- cbind(X,"Age" = (covar.df$Age - mean(covar.df$Age)) / sd(covar.df$Age))
    }

    if ("sex" %in% covariates) {
      X <- cbind(X,"sex:Male vs. Female" = 1 * (covar.df$sex == "M"))
    }

    if ("eth" %in% covariates) {
      X <- cbind(X,"eth:afr vs. eur" = 1 * (covar.df$eth == "afr"), "eth:asi vs. eur" = 1 * (covar.df$eth == "asi"))
    }
    if ("treatment" %in% covariates) {
      X <- cbind(X,"treatment:Test vs. Placebo" = 1 * (covar.df$treatment == "Test"))
    }
    if ("severity" %in% covariates) {
      X <- cbind(
        X,
        "severity:Low vs. Very Low" = 1 * (covar.df$severity == "Low"),
        "severity:Moderate vs. Very Low" = 1 * (covar.df$severity == "Moderate"),
        "severity:High vs. Very Low" = 1 * (covar.df$severity == "High"),
        "severity:Very High vs. Very Low" = 1 * (covar.df$severity == "Very High"),
        "severity:Extreme vs. Very Low" = 1 * (covar.df$severity == "Extreme")
      )
    }

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
  
  names(pheno.df)[names(pheno.df) == "ind"] <- "IID"
  pheno.plink.df <- data.frame("#FID" = fid.id, pheno.df, check.names = FALSE)

  write.table(pheno.plink.df[, c(1, 2, 3)], file = num.pheno.plink.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  write.table(pheno.plink.df[, c(1, 2, 4)], file = bin.pheno.plink.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
  

  remove(list = "pheno.df")
  if (verbose) cat("generated phenotype files\n")
  num.poi.blocks <- ceiling(num.poi / poi.chunk.size)
  if (poi.file.type == "h5") {
    dir.create(poi.data.dir, showWarnings = FALSE)
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
    
    # maf <- runif(num.poi, min = 0.05, max = 0.5)
    # miss.rate <- runif(num.poi, min = 0, max = 0.1)
    # # for (chunk_start in seq(1, num.poi, by=poi.chunk.size)) {
    #   # chunk_end <- min(chunk_start + poi.chunk.size - 1, num.poi)
    # chunk_indices <- 1:num.poi
    # values <- generate_values(num.ind, chunk_indices, poi.type, maf[chunk_indices], miss.rate[chunk_indices])
    # colnames(values) <- poi.id[chunk_indices]
    # df <- data.frame(ID = ind.id, values, stringsAsFactors = FALSE, check.names = FALSE)
    # write.table(df, file=poi.txt.file, sep="\t", row.names = FALSE, col.names = (chunk_indices[1] == 1), append = (chunk_indices[1] != 1), na = "", quote = FALSE)
    # # fwrite(df, file = poi.txt.file, sep = "\t", row.names = FALSE, col.names = (chunk_indices[1] == 1), append = (chunk_indices[1] != 1), na = "", quote = FALSE)
    # # }
    # # values <- generate_values(num.ind, 1:num.poi, poi.type, poi.chunk.size)
    # # colnames(values) <- poi.id
    # # values <- data.frame(ID=ind.id, values, stringsAsFactors = FALSE, check.names = FALSE)
    # # fwrite(values, file = poi.txt.file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
    # if (verbose) cat("generated poi txt file\n")
    # Step 1: Write the ID column to a file
write.table(data.frame(ID = ind.id), file = "ID_column.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Step 2: Set up parallel processing
num_cores <- detectCores() - 4  # Use four less than the total cores to prevent overloading
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Step 3: Calculate the total number of chunks
num_chunks <- ceiling(num.poi / poi.chunk.size)
chunk_nums <- 1:num_chunks

# Step 4: Generate overall 'maf' and 'miss.rate' vectors
maf <- runif(num.poi, min = 0.05, max = 0.5)
miss.rate <- runif(num.poi, min = 0, max = 0.1)

# Step 5: Process chunks in parallel
data_files_list <- foreach(chunk_num = chunk_nums, .packages = c()) %dopar% {
  # Calculate chunk indices
  chunk_start <- (chunk_num - 1) * poi.chunk.size + 1
  chunk_end <- min(chunk_num * poi.chunk.size, num.poi)
  chunk_indices <- chunk_start:chunk_end
  num.poi_chunk <- length(chunk_indices)

  # Extract 'maf' and 'miss.rate' for the current chunk
  maf_chunk <- maf[chunk_indices]
  miss_rate_chunk <- miss.rate[chunk_indices]

  # Generate values
  values <- generate_values(num.ind, num.poi_chunk, poi.type, maf_chunk, miss_rate_chunk)

  # Set column names
  colnames(values) <- poi.id[chunk_indices]

  # Write only the data columns to a temporary file
  temp_file <- paste0("chunk_data_", chunk_num, "_", Sys.getpid(), ".txt")
  write.table(values, file = temp_file, sep = "\t", row.names = FALSE, col.names = TRUE,
              na = "", quote = FALSE)

  if (verbose) cat("Processed chunk", chunk_num, "of", num_chunks, "\n")

  return(temp_file)
}

# Step 6: Stop the cluster after processing
stopCluster(cl)

# Step 7: Unlist the results to get a character vector of file names
data_files <- unlist(data_files_list)

# Include the ID column file in the list of data files
data_files <- c("ID_column.txt", data_files)

# Step 8: Prepare to merge data files
# Open the output file for writing
output_conn <- file(poi.txt.file, open = "w")

# Write the header line
header_line <- c("ID", poi.id)
writeLines(paste(header_line, collapse = "\t"), output_conn)

rows_per_chunk <- 1000  # Adjust based on memory capacity

# Initialize a progress counter
if (verbose) cat("Merging data files...\n")

# Initialize positions for data files (excluding ID column)
data_positions <- vector("list", length = length(data_files) - 1)

# Open the ID column file
id_conn <- file("ID_column.txt", open = "r")
# Read and discard the header line
readLines(id_conn, n = 1)

# Initialize line counter
line_num <- 0

# Step 9: Merge data files incrementally
repeat {
  # Read a chunk of IDs
  id_lines <- readLines(id_conn, n = rows_per_chunk)

  # Break if no more lines
  if (length(id_lines) == 0) break

  # Convert ID lines to a data frame
  id_chunk <- read.table(text = id_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

  # Initialize a list to hold data chunks
  data_chunks <- list()

  # Process each data file (excluding ID column)
  for (i in seq_along(data_files)[-1]) {
    temp_file <- data_files[[i]]  # Use double brackets to extract the filename

    # Open the data file if not already open
    if (is.null(data_positions[[i - 1]]$conn)) {
      data_positions[[i - 1]] <- list(
        conn = file(temp_file, open = "r"),
        eof = FALSE
      )
      # Read and discard the header line
      readLines(data_positions[[i - 1]]$conn, n = 1)
    }

    # Read a chunk of data lines
    data_lines <- readLines(data_positions[[i - 1]]$conn, n = rows_per_chunk)

    # Check for end of file
    if (length(data_lines) == 0) {
      data_positions[[i - 1]]$eof <- TRUE
      # Close the connection
      close(data_positions[[i - 1]]$conn)
      # Remove the temporary file to save space
      file.remove(temp_file)
      next
    }

    # Convert data lines to a data frame
    data_chunk <- read.table(text = data_lines, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

    # Add to the list of data chunks
    data_chunks[[length(data_chunks) + 1]] <- data_chunk
  }

  # Combine the ID chunk and data chunks
  combined_chunk <- cbind(id_chunk, do.call(cbind, data_chunks))

  # Write the combined chunk to the output file
  write.table(combined_chunk, file = output_conn, sep = "\t", row.names = FALSE,
              col.names = FALSE, na = "", quote = FALSE, append = TRUE)

  # Update line counter
  line_num <- line_num + nrow(combined_chunk)

  if (verbose) cat("Merged", line_num, "rows\n")
}

# Step 10: Close file connections and clean up
# Close the ID connection
close(id_conn)
# Close the output connection
close(output_conn)
# Remove the ID column file
file.remove("ID_column.txt")

if (verbose) cat("All chunks processed and merged into", poi.txt.file, "\n")
  }

  # file.con[["values"]] <- t(values);
  # file.con$close_all();

  invisible(list("num.poi" = num.poi, "num.ind" = num.ind, "num.resp.beta0" = num.resp.beta0, "num_resp.beta" = num_resp.beta, "num.resp.sd" = num.resp.sd, "bin.resp.beta0" = bin.resp.beta0, "bin_resp.beta" = bin_resp.beta))
}

#' generate_values function to generate values for the test dataset
#' @param num.ind an integer denoting number of individuals
#' @param poi.type character must be either "genotype" or "dosage"
#' @param maf MAF threshold 
#' @param miss.rate rate for missing values 
#' @return a list consisting of dataset size 
#' @import stats
#' @import parallel
#' @import data.table
#' @import rhdf5
#' @import utils
# Define the generate_values function
generate_values <- function(num.ind, num.poi, poi.type, maf, miss.rate) {
  if (poi.type == "genotype") {
    # Generate dosage values and thresholds
    dosage.val <- matrix(runif(num.ind * num.poi), nrow = num.ind, ncol = num.poi)
    maf.mat <- matrix(maf, nrow = num.ind, ncol = num.poi, byrow = TRUE)
    t0.mat <- (1 - maf.mat) ^ 2
    t1.mat <- 1 - maf.mat ^ 2

    # Initialize geno.val
    geno.val <- matrix(NA_real_, nrow = num.ind, ncol = num.poi)
    
    # Assign genotype values based on dosage thresholds
    geno.val[dosage.val < t0.mat] <- 0
    geno.val[dosage.val >= t0.mat & dosage.val < t1.mat] <- 1
    geno.val[dosage.val >= t1.mat] <- 2
    
    # Introduce missingness
    miss.val <- matrix(runif(num.ind * num.poi), nrow = num.ind, ncol = num.poi)
    miss.rate.mat <- matrix(miss.rate, nrow = num.ind, ncol = num.poi, byrow = TRUE)
    geno.val[miss.val > 1 - miss.rate.mat] <- NA
    
    values <- geno.val
  } else {
    # Generate random values for non-genotype data
    r0 <- matrix(runif(num.ind * num.poi, min = 0, max = 0.4), nrow = num.ind)
    r1 <- matrix(runif(num.ind * num.poi, min = 0.25, max = 0.75), nrow = num.ind)
    r2 <- matrix(runif(num.ind * num.poi, min = 0.6, max = 1.0), nrow = num.ind)
    
    maf.mat <- matrix(maf, nrow = num.ind, ncol = num.poi, byrow = TRUE)
    c0.mat <- (1 - maf.mat) ^ 2
    c1.mat <- 2 * (1 - maf.mat) * maf.mat
    c2.mat <- maf.mat ^ 2
    
    # Compute genotype values
    geno.val <- c0.mat * r0 + c1.mat * r1 + c2.mat * r2
    
    # Introduce missingness
    miss.val <- matrix(runif(num.ind * num.poi), nrow = num.ind)
    miss.rate.mat <- matrix(miss.rate, nrow = num.ind, ncol = num.poi, byrow = TRUE)
    geno.val[miss.val > 1 - miss.rate.mat] <- NA
    
    values <- geno.val
  }
  
  values
}


#' generate_values2 function to generate values for the test dataset
#' @param num.ind an integer denoting number of individuals
#' @param poi.indices vector of POI ids
#' @param poi.type character must be either "genotype" or "dosage"
#' @param poi.chunk.size chunk size for the hdf5 dataset
#' @return a list consisting of dataset size 
#' @import stats
#' @import parallel
#' @import data.table
#' @import rhdf5
#' @import utils
generate_values2 <- function(num.ind, poi.indices, poi.type, poi.chunk.size) {
  block.index <- poi.indices
  block.size <- length(block.index)
  # cat("block.size: ")
  # cat(block.size)
  # cat("num.ind: ")
  # cat(num.ind)
  # miss.rate <- runif(block.size, min = 0, max = 0.1)

  if (poi.type == "genotype") {
    values <- matrix(0.000, ncol = block.size, nrow = num.ind)
  } else {
    values <- matrix(0.000, ncol = block.size, nrow = num.ind)
  }
  maf <- runif(block.size, min = 0.05, max = 0.5)
  miss.rate <- runif(block.size, min = 0.0, max = 0.1)
  for (i in 1:block.size) {
    if (poi.type == "genotype") {
      dosage.val <- runif(num.ind)
      geno.val <- integer(num.ind)
      geno.val[dosage.val < (1.0 - maf[i])^2] <- 0.0
      geno.val[((1 - maf[i])^2 < dosage.val) & (dosage.val < (1 - maf[i]^2))] <- 1.0
      geno.val[dosage.val > (1.0 - maf[i]^2)] <- 2.0
      geno.val[runif(num.ind)> 1.0 - miss.rate[i]] <- NA;
      values[, i] <- geno.val
    } else {
      geno.val <- ((1 - maf[i])^2) * runif(num.ind, min = 0, max = 0.4) + 2.0 * (1.0 - maf[i]) * maf[i] * runif(num.ind, min = 0.25, max = 0.75) + (maf[i]^2) * runif(num.ind, min = 0.6, max = 1.0)
      geno.val[runif(num.ind) > 1 - miss.rate[i]] <- NA
      values[, i] <- geno.val
      # miss.rate <- runif(block.size, min = 0, max = 0.1)
    }
  }
  
  values
}

#' write.h5 function to write the hdf5 dataset
#' @param data.num.chunks number of files to split the data into
#' @param prefix prefix for all files created
#' @param verbose logical to control display of progress messages (default=TRUE)
#' @param ind.id array of unique individual ids
#' @param poi.id array of unique POI ids
#' @param poi.data.dir output directory for hdf5 dataset
#' @param poi.compression.level gzip compression level for hdf5
#' @param poi.chunk.size chunk size for the hdf5 dataset
#' @param poi.type character must be either "genotype" or "dosage"
#' @param num.ind an integer denoting number of individuals
#' @param num.poi an integer denoting number of predictors of interest (SNP calls)
#' @return none 
#' @import stats
#' @import parallel
#' @import data.table
#' @import rhdf5
#' @import utils
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
      storage.mode = "double",
      chunk = c(num.ind, min(poi.chunk.size, pois.per.file[k])),
      level = poi.compression.level
    )
    for (block.start in seq(start.index, end.index, by = poi.chunk.size)) {
      block.end <- min(block.start + poi.chunk.size - 1, end.index)
      values <- generate_values2(num.ind, block.start:block.end, poi.type)
      h5write(values, file = poi.file, name = "values", start = c(1, block.start - start.index + 1), count = c(num.ind, block.end - block.start + 1))
    }
  }
  if (verbose) cat("added POI values into h5 file\n")
  h5closeAll()
}
