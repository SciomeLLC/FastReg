library("data.table");
library("hdf5r");
library("FastReg");

setwd("C:/Users/deepak.mav/Documents/data/Projects/NIEHS_Biostat/Alison-Motsinger/FastReg/Packages/data");

num.out <- FastReg(config.file="C:/Users/deepak.mav/Documents/data/Projects/NIEHS_Biostat/Alison-Motsinger/FastReg/Packages/data/example.num.config")
#bin.out <- FastReg(config.file="example.bin.config")
