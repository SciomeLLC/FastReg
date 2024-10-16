setwd("../fastvla_data")
library(FastVLA)
#aa = FastVLA::generate.dataset(seed = 123, num.pc = 20, num.subj = 500, num.variants = 5000, num.pheno = 5)
testt = BEDMatrix::BEDMatrix("plink/testdata_5k_by_50k.bed")

Y = readRDS("Y_log.rds")
X = readRDS("X.rds")
a = proc.time()
FastVLA_logistic(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_single_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 1, max_blas_threads = 1)
show(proc.time() - a)
a = proc.time()
FastVLA_logisticf(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_float_single_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 1, max_blas_threads = 1)
show(proc.time() - a)

# Multi threaded
a = proc.time()
FastVLA_logistic(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_two_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 2, max_blas_threads = 1)
show(proc.time() - a)
a = proc.time()
FastVLA_logisticf(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_float_two_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 2, max_blas_threads = 1)
show(proc.time() - a)

# Multi threaded, multi blas threaded
# a = proc.time()
# FastVLA_logistic(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_multi_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 1, max_blas_threads = 4)
# show(proc.time() - a)
# a = proc.time()
# FastVLA_logisticf(Y, slot(testt, "xptr"), 1:5000, 1:500, X, 100, "test_float_multi_thread", colnames(Y), paste0("Y", 1:5000), "VB_MR2", epss = 1e-16, mafthresh = 0.08, max_iter = 10, max_threads = 1, max_blas_threads = 4)
# show(proc.time() - a)