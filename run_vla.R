library(FastReg)
setwd("../jsh.input/Input/fastvla")

mats <- VLA(
  regression.type="logistic", # logistic or linear
  pheno.file="testdata_100_by_1000.bin.pheno.txt",  # phenotype file to use
  pheno.rowname.cols="ind", # ID column that matches rows of the HDF5 file(s) created
  phenotype = "bin.resp", # column name of phenotype (trait) in pheno file
  covar.file="testdata_100_by_1000.covar.txt", # covariate file to use
  covar.rowname.cols="ind", # ID column that matches rows of the HDF5 file(s) created
  covariates=c("age", "sex", "treatment"),
  covariate.type = c("numeric", "categorical", "categorical"),
  covariate.standardize = c(FALSE, FALSE, FALSE),
  covariate.levels = c("", "M,F", "Test,Placebo"),
  covariate.ref.level = c("", "F", "Test"),
  POI.file.dir = ".",
  POI.type = "genotype", # unless integer, keep as dosage
  POI.file.format = "bed",
  POI.effect.type = "additive", # keep as additive unless you have integer SNP data and are calling a dominant or recessive model,
  POI.covar.interactions = c("age", "treatment"),
  output.dir = "win_R4.2.2")