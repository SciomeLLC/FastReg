source("simulate_test_data_function.R");

simulate_test_dataset(num.poi = 50000,
                                  num.ind = 5000,
                                  seed = 12133,
                                  coeff.sd = 0,
                                  bin.resp.mean = 0.2,
                                  num.resp.mean = 24,
                                  num.resp.sd = 5,
                                  poi.type = "genotypes",
                                  poi.chunk.size = 1000,
                                  poi.compression.level=7,
                                  data.dir = "c:/users/deepak.mav/Documents/data/temp",
                                  prefix = "testdata_5k_by_50k", verbose=TRUE)