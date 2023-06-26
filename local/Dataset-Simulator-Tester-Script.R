source("simulate_test_data_function.R");

simulate_test_dataset(num.poi = 100000,
                                  num.ind = 5000,
                                  seed = 12133,
                                  coeff.sd = 0,
                                  bin.resp.mean = 0.1,
                                  num.resp.mean = 10,
                                  num.resp.sd = 5,
                                  poi.type = "genotypes",
                                  poi.chunk.size = 10000,
                                  poi.compression.level=7,
                                  data.dir = "F:/",
                                  prefix = "testdata_10p_5k_x_10m", verbose=TRUE)
