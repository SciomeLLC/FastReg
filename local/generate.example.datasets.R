num.poi <- 50;
num.ind <- 5;
covariates <- c("age", "sex", "eth", "treatment", "severity");

seed <- 12133;
poi.type <- "genotypes";
poi.chunks <- 100;
data.dir <- "E:/Development";
setwd(data.dir);

covar.file <- "small.covar.txt";
num.pheno.file <- "small.num.pheno.txt";
bin.pheno.file <- "small.bin.pheno.txt";
poi.file <- "small.poi.h5";
poi.subset.file <- "small.poi.subset.txt"
subject.subset.file <- "small.sample.subset.txt"


set.seed(seed);

ind.id <- paste0("IND", formatC(1:num.ind, format="d", flag="0", digits=floor(log10(num.ind))));
poi.id <- unique(paste0("rs", formatC(10000000*runif(2*num.poi), format="d", flag="0", digits=8)))[1:num.poi];


write.table(data.frame("ind"=sort(sample(ind.id, size=5, replace=FALSE)), stringsAsFactors=FALSE), file=subject.subset.file, quote=FALSE, na="", row.names=FALSE, col.names=TRUE);
write.table(data.frame("rsID"=sort(sample(poi.id, size=50, replace=FALSE)), stringsAsFactors=FALSE), file=poi.subset.file, quote=FALSE, na="", row.names=FALSE, col.names=TRUE);


covar.df <- data.frame("ind"=ind.id,  "age" = as.integer(rnorm(num.ind, 50, 7)),
                          "sex"=c("M","F")[ceiling(2*runif(num.ind))],
                          "eth" = c("eur", "afr", "asi")[ceiling(3*runif(num.ind))],
                          "treatment" = c("Test", "Placebo")[ceiling(2*runif(num.ind))],
                          "severity" = c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme")[ceiling(6*runif(num.ind))],
                          stringsAsFactors=FALSE, check.names=FALSE);

write.table(covar.df, file=covar.file, sep="\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
remove(list="covar.df");

pheno.df <- data.frame("ind"=ind.id, "num.resp"=rnorm(num.ind,24, 5), "bin.resp"=1*(runif(num.ind)>=0.8));

write.table(pheno.df[, c(1,2)], file=num.pheno.file, sep="\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
write.table(pheno.df[, c(1,3)], file=bin.pheno.file, sep="\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
remove(list="pheno.df");


library("hdf5r");
file.con <- H5File$new(poi.file, mode="w");
file.con[["individuals"]] <- ind.id;
file.con[["predictors_of_interest"]] <- poi.id;

maf <- runif(num.poi, min=0.05, max=0.5);
miss.rate <- runif(num.poi, min=0, max=0.1);

values <- matrix(0, ncol=num.poi, nrow=num.ind);
for(i in 1:num.poi) {
  dosage.val <- runif(num.ind);
  geno.val <- integer(num.ind);
  geno.val[dosage.val < (1-maf[i])^2] <- 0;
  geno.val[((1-maf[i])^2 < dosage.val) & (dosage.val < (1 - maf[i]^2))] <- 1;
  geno.val[dosage.val > (1- maf[i]^2)] <- 2;
  geno.val[runif(num.ind)>1-miss.rate[i]] <- NA;
  values[,i] <- geno.val
}

file.con[["values"]] <- t(values);
file.con$close_all();
remove(list=c("dosage.val","geno.val", "maf", "miss.rate"))



