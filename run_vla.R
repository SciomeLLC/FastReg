setwd("../fastvla_data")
library(FastReg)
#aa = FastVLA::generate.dataset(seed = 123, num.pc = 20, num.subj = 500, num.variants = 5000, num.pheno = 5)
testt = BEDMatrix::BEDMatrix("plink/testdata_5k_by_50k.bed")

Y = readRDS("Y_log.rds")

X = readRDS("X.rds")
a = proc.time()
FastVLA_chunked_sota(Y, slot(testt, "xptr"), 1:50, 1:500, X, 100, "test2", colnames(Y), paste0("Y", 1:100), "VB_MR4", epss = 1e-16)
show(proc.time() - a)