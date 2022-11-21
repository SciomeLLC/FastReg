library("rhdf5") # bioconductor
library(parallel)
library(memuse)

# pass arguments for max_cores to use
get_chunk_size <- function(hdf_file_name) {
  has_blas <- TRUE
  cores <- detectCores(logical = TRUE)
  os <- Sys.info()[["sysname"]]
  # MacOS
  if (os == "Darwin") {
    memfree <- as.numeric(system("top -l1 -s0 | awk '/PhysMem/ {print $6+0}'", intern = TRUE)) * 1024 * 1024 # convert to bytes
  } else {
    memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024 # convert to bytes
  }

  # check for BLAS
  #if (flexiblas_avail()) {
  #  blas_libs <- flexiblas_list()
  #  has_blas <- length(blas_libs) > 0
  #}

  float_size <- 8 # 8 bytes per number assuming 64-bit numbers
  listFile <- h5ls(hdf_file_name)
  dim1 <- as.numeric(listFile$"dim"[[1]])
  dim2 <- as.numeric(listFile$"dim"[[2]])

  data_size <- dim1 * dim2 * float_size
  # shave off 2 workers / 1 core
  chunks <- data_size / (memfree * 0.8)
  chunked_dim1 <- floor(dim1 / chunks)
  chunk_size <- list(chunked_dim1, dim2)
  chunked_parallel <- floor(chunked_dim1 / cores)

  if (has_blas == TRUE) {
    return(list(chunked_parallel, cores, memfree))
  } else {
    return(list(chunked_dim1, 1, memfree))
  }
}

chunk_size <- get_chunk_size("random.h5")

available_cores <- chunk_size[[2]]
chunks <- chunk_size[[1]]
memfree <- chunk_size[[3]]

print(memfree)
dat <- h5read("random.h5", "/values", index = list(NULL, 1:chunks))

dat
dat <- h5read("random.h5", "/values", index = list(NULL, chunks + 1:chunks * 2))
remove(dat)
h5closeAll()
gc()
