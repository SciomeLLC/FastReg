source("simulate_test_data_function.R");

simulate_test_dataset(num.poi = 5000,
                                  num.ind = 5000,
                                  seed = 12133,
                                  coeff.sd = 0,
                                  bin.resp.mean = 0.2,
                                  num.resp.mean = 24,
                                  num.resp.sd = 5,
                                  poi.type = "additive",
                                  poi.chunk.size = 1000,
                                  poi.compression.level=7,
                                  data.dir = "F:/",
                                  prefix = "testdata_5000_x_5000_new", verbose=TRUE)