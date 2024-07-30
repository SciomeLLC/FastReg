## ----include=TRUE-------------------------------------------------------------
# Binary phenotype data
read.csv("../data/testdata.bin.pheno.txt", sep = "\t", nrows = 10)

## ----include=TRUE-------------------------------------------------------------
# Covariate data
read.csv("../data/testdata.covar.txt", sep = "\t", nrows = 10)

## ----eval=TRUE, include=TRUE--------------------------------------------------
# POI data
read.csv("../data/testdata.poi.txt", sep = "\t", nrows = 2)

## ----include=TRUE, eval=FALSE-------------------------------------------------
#  library(FastReg)
#  FastRegImport("testdata", "./FastRegData", delimiter="\t")

## ----include=TRUE, eval=FALSE-------------------------------------------------
#  FastReg(
#    regression.type="logistic",
#    pheno.file="data/testdata.bin.pheno.txt",
#    pheno.rowname.cols="ind",
#    phenotype = "bin.resp",
#    covar.file="data/testdata.covar.txt",
#    covariates=c("age", "sex", "treatment"),
#    covariate.type = c("numeric", "categorical", "categorical"),
#    covariate.standardize = c(TRUE, FALSE, FALSE),
#    covariate.levels = c("", "F,M", "Placebo,Test"),
#    covariate.ref.level = c("", "F", "Placebo"),
#    POI.file.dir = "./FastRegData"
#    POI.type = "dosage",
#    POI.effect.type = "additive",
#    output.dir = "FastRegResults"
#  )

## ----include=TRUE, eval=FALSE-------------------------------------------------
#  read.csv("FastRegResults/Full_Convergence_stratum_1.tsv", delimiter="\t", nrows=10)
#  read.csv("FastRegResults/Full_Results_stratum_1.tsv", delimited="\t", nrwos=10)

