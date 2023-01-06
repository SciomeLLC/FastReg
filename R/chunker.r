estimate.poi.block.size <- function(num.poi, num.ind, poi.type, max.cores = NULL) {
  num_threads <- detectCores(logical = TRUE)
  if (!is.null(max.cores) && num_threads > max.cores) {
    num_threads <- max.cores
  }
  os <- Sys.info()[["sysname"]]

  if (os == "Windows") {
    x <- system2("wmic", args = "OS get FreePhysicalMemory /Value", stdout = TRUE)
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    memfree <- as.numeric(x) * 1024 # convert to bytes
  } else if (os == "Darwin") { # MacOS
    memfree <- as.numeric(system("top -l1 -s0 | awk '/PhysMem/ {print $6+0}'", intern = TRUE)) * 1024 * 1024 # convert to bytes
  } else { # Linux
    memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024 # convert to bytes
  }

  # keep one thread idle
  if (num_threads > 1) {
    num_threads <- num_threads - 1
  }

  matrix_size <- exp(log(num.poi) + log(num.ind))
  float_size <- 8L # 8 bytes per number assuming 64-bit numbers
  data_size <- exp(log(matrix_size) + log(as.numeric(float_size)) + log(40L))
  master_thread_memory <- 524288000L # 500mb

  chunks <- (data_size + master_thread_memory) / ((memfree) * 0.8)
  chunked_dim1 <- floor(num.poi / chunks)
  # chunk_size <- list(chunked_dim1, num.ind)

  chunked_parallel <- floor(chunked_dim1 / num_threads)
  # stop("No overlapping individuals found in POI, pheno, covar files");
  return(chunked_parallel)
}
