library("hdf5r");
#library("flexiblas");
library("data.table");


setwd("C:/Users/deepak.mav/Documents/data/Projects/NIEHS_Biostat/Alison-Motsinger/FastReg/Packages/FastReg/R");
source("Array_Functions.R");
source("logistic.R");
source("linear.R")
source("util-functions.R");
source("parserFunctions.R");
source("main.R")

setwd("C:/Users/deepak.mav/Documents/data/Projects/NIEHS_Biostat/Alison-Motsinger/FastReg/Packages/FastReg/data");

num.out <- FastReg(config.file="example.num.config")
#bin.out <- FastReg(config.file="example.bin.config")
