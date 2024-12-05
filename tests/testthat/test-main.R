test_that("FastVLA_single_Y works as expected", {
  library(BEDMatrix)

  # Load the test data
  test_data_dir <- testthat::test_path("testdata")
  X <- readRDS(file.path(test_data_dir, "newX.rds"))
  Y <- readRDS(file.path(test_data_dir, "newY.rds"))
  G_pl <- BEDMatrix::BEDMatrix(file.path(test_data_dir, "plinkfile.bed"))

  # Run the function and measure execution time
  start_time <- proc.time()
  result <- FastVLA::FastVLA_single_Y(
    c(Y[, 1]),
    slot(G_pl, "xptr"),
    1:5000,
    1:5000,
    X,
    10,
    ".",
    c("Y1", "Y2"),
    paste0("POI", 1:5000),
    "dub",
    epss = 1e-16,
    mafthresh = 0.01,
    do_pca = FALSE,
    add_intercept = FALSE
  )
  end_time <- proc.time() - start_time
  # Show the execution time
  message("Execution time: ", end_time)
})