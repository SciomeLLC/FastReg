source("simulate_test_data_function.R");

simulate_test_dataset(num.poi = 500,
                                  num.ind = 500,
                                  seed = 12133,
                                  coeff.sd = 0,
                                  bin.resp.mean = 0.2,
                                  num.resp.mean = 24,
                                  num.resp.sd = 5,
                                  poi.type = "genotypes",
                                  poi.chunk.size = 1000,
                                  poi.compression.level=7,
                                  data.dir = "F:/",
                                  prefix = "testdata_500_x_500", verbose=TRUE)
