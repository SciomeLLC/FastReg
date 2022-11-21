library("hdf5r");
# library("flexiblas");
library("data.table");


setwd("E:/Development/FastR/FastReg/R");
source("Array_Functions.R");
source("logistic.R");
source("linear.R")
source("util-functions.R");
source("parserFunctions.R");
source("main.R")


setwd("E:/Development/FastR/FastReg/data");

# bin.out <- FastReg(config.file="example.bin.config")
num.out <- FastReg(config.file="example.num.config")

