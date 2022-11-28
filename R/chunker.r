library("rhdf5") # bioconductor
library(parallel)
library(memuse)
library(flexiblas)

# pass arguments for max_cores to use
estimate.poi.block.size <- function(num.poi, num.ind, poi.type, num.cores) {
  has_blas <- TRUE
  cores <- detectCores(logical = TRUE)
  os <- Sys.info()[["sysname"]]
  # MacOS
  if (os == "Darwin") {
    memfree <- as.numeric(system("top -l1 -s0 | awk '/PhysMem/ {print $6+0}'", intern = TRUE)) * 1024 * 1024 # convert to bytes
  } else {
    memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE)) * 1024 # convert to bytes
  }

  # check for BLAS and multi-threading
  if (flexiblas_avail()) {
    blas_libs <- flexiblas_list()
    has_blas <- length(blas_libs) > 0
    num_threads <- flexiblas_get_num_threads()
    idx <- flexiblas_load_backend(backends)

    blas_lib_id <- 1

    for (id in idx){
      flexiblas_switch(id)
      num_threads_temp <- flexiblas_get_num_threads()
      if (num_threads_temp > num_threads){
        num_threads = num_threads_temp
        blas_lib_id <- id
      }
    }
  }
  # switch to blas library with highest threads capability
  flexiblas_switch(blas_lib_id)
  if (!is.null(num.cores) && num.cores *2 < num_threads) {
    num_threads = num.cores * 2
  }

  float_size <- 8 # 8 bytes per number assuming 64-bit numbers
  data_size <- num.poi * num.ind * float_size
  # shave off 2 workers / 1 core
  chunks <- data_size / (memfree * 0.8)
  chunked_dim1 <- floor(num.poi / chunks)
  chunk_size <- list(chunked_dim1, num.ind)
  chunked_parallel <- floor(chunked_dim1 / num_threads)

  return(chunked_parallel)
}

# chunk_size <- get_chunk_size("random.h5")

# available_cores <- chunk_size[[2]]
# chunks <- chunk_size[[1]]
# memfree <- chunk_size[[3]]

# print(memfree)
# dat <- h5read("random.h5", "/values", index = list(NULL, 1:chunks))

# dat
# dat <- h5read("random.h5", "/values", index = list(NULL, chunks + 1:chunks * 2))
# remove(dat)
# h5closeAll()
# gc()
