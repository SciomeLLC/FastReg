library("rhdf5")
library("parallel")
#' simulate_test_dataset function to generate test dataset for evaluation of FastReg
#'@param num.poi an integer denoting number of predictors of interest (SNP calls)
#'@param num.ind an integer denoting number of individuals
#'@param seed an integer denoting random number generate seed
#'@param coeff.sd numeric value denoting the standard deviation parameter for regression coefficients for covariates (default=0, regression coefficients for covarites are set to zero)
#'@param bin.resp.mean numeric value between 0 and 1 denoting baseline incidence rate for binary response
#'@param num.resp.mean numeric value denoting overall mean of numeric response
#'@param num.resp.sd numeric value denoting overall standard deviation of numeric response
#'@param poi.type character must be either "genotypes" or "dosage"
#'@param poi.chunk.size an integer denoting poi chunk size used during H5 file generation
#'@param poi.compression.level an integer denoting compression level used during H5 file creation
#'@param verbose logical to control display of progress messages (default=TRUE)
#'@return a list consisting of dataset size (num.poi, num.ind) as well as regression coefficients used for both binary and numeric response
simulate_test_dataset <- function(num.poi = 50000,
                                  num.ind = 5000,
                                  seed = 12133,
								                  coeff.sd = 0,
								                  bin.resp.mean = 0.2,
								                  num.resp.mean = 24,
								                  num.resp.sd = 5,
								                  poi.type = "genotype",
								                  poi.chunk.size = 100,
								                  poi.compression.level = 7,
                                  data.dir = ".",
                                  prefix = "testdata_5k_by_50k", 
                                  verbose=TRUE,
                                  poi.file.type = "txt"){
  # num.poi <- 50000
  # num.ind <- 5000
  covariates <- c("age", "sex", "eth", "treatment", "severity")

  # seed <- 12133;
  #poi.type <- "genotypes"
  #poi.chunk.size <- 100
  #poi.compression.level <- 7
  # data.dir <- "../input/"
  # prefix = "testdata_5k_by_50k"

  dir.create(data.dir, showWarnings=FALSE);

  covar.file <-     file.path(data.dir, paste0(prefix, ".covar.txt"));
  num.pheno.file <- file.path(data.dir, paste0(prefix, ".num.pheno.txt"));
  bin.pheno.file <- file.path(data.dir, paste0(prefix, ".bin.pheno.txt"));
  poi.txt.file <-       file.path(data.dir, paste0(prefix, ".poi.txt"));
  poi.file <-       file.path(data.dir, paste0(prefix, ".poi.h5"));
  poi.subset.file <- file.path(data.dir, paste0(prefix, ".poi.subset.txt"));
  subject.subset.file <-file.path(data.dir, paste0(prefix, ".sample.subset.txt"));


  set.seed(seed);

  ind.id <- paste0("IND", formatC(1:num.ind, format="d", flag="0", digits=floor(log10(num.ind))));
  poi.id <- paste0("rs", formatC(1:num.poi, format="d", flag="0", digits=floor(log10(num.poi))));


  write.table(data.frame("ind"=sort(sample(ind.id, size=5, replace=FALSE)), stringsAsFactors=FALSE),
              file=subject.subset.file, quote=FALSE, na="", row.names=FALSE, col.names=TRUE);
  write.table(data.frame("rsID"=sort(sample(poi.id, size=50, replace=FALSE)), stringsAsFactors=FALSE),
              file=poi.subset.file, quote=FALSE, na="", row.names=FALSE, col.names=TRUE);


  covar.df <- data.frame("ind"=ind.id,  "age" = as.integer(rnorm(num.ind, 50, 7)),
                         "sex"=c("M","F")[ceiling(2*runif(num.ind))],
                         "eth" = c("eur", "afr", "asi")[ceiling(3*runif(num.ind))],
                         "treatment" = c("Test", "Placebo")[ceiling(2*runif(num.ind))],
                         "severity" = c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme")[ceiling(6*runif(num.ind))],
                         stringsAsFactors=FALSE, check.names=FALSE);

  write.table(covar.df, file=covar.file, sep="\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
  remove(list="covar.df");
  if(verbose) cat("generated covariate file\n");

  num.resp.beta0 <- num.resp.mean;
  bin.resp.beta0 <-  log(bin.resp.mean/(1-bin.resp.mean));


  if(coeff.sd>0) {
	X <- cbind("Age"=(covar.df$Age-mean(covar.df$Age))/sd(covar.df$Age),
					"sex:Male vs. Female"=1*(covar.df$sex=="M"),
					"eth:afr vs. eur"= 1*(covar.df$eth=="afr"), "eth:asi vs. eur"=1*(covar.df$eth=="asi"),
					"treatment:Test vs. Placebo"=1*(covar.df$treatment=="Test"),
					"severity:Low vs. Very Low"=1*(covar.df$severity=="Low"),
					"severity:Moderate vs. Very Low"=1*(covar.df$severity=="Moderate"),
					"severity:High vs. Very Low"=1*(covar.df$severity=="High"),
					"severity:Very High vs. Very Low"=1*(covar.df$severity=="Very High"),
					"severity:Extreme vs. Very Low"=1*(covar.df$severity=="Extreme"));


	num_resp.beta <- matrix(rnorm(ncol(X), mean=0, sd=coeff.sd), ncol=1, dimnames=list(colnames(X),NULL));
	bin_resp.beta <- matrix(rnorm(ncol(X), mean=0, sd=coeff.sd), ncol=1, dimnames=list(colnames(X),NULL));
	num.resp <- rnorm(num.ind, mean=num.resp.beta0 + X%*%num_resp.beta, sd=num.resp.sd);
	bin.resp <- rbinom(num.ind, prob=1/(1+exp(-(bin.resp.beta0+X%*%bin_resp.beta))),size=1);
	remove(list="X");
 } else {
	num.resp <- rnorm(num.ind, mean=num.resp.beta0, sd=num.resp.sd);
	bin.resp <- rbinom(num.ind, prob=1/(1+exp(-(bin.resp.beta0))), size=1);
	num_resp.beta <- NULL;
	bin_resp.beta <- NULL;
 }

  pheno.df <- data.frame("ind"=ind.id, "num.resp"=num.resp, "bin.resp"=bin.resp);


  write.table(pheno.df[, c(1,2)], file=num.pheno.file, sep = "\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
  write.table(pheno.df[, c(1,3)], file=bin.pheno.file, sep = "\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE);
  remove(list="pheno.df");
  if(verbose) cat("generated phonotype files\n");
  if(poi.file.type == "h5") {
    if(file.exists(poi.file)) unlink(poi.file);
    #file.con <- H5File$new(poi.file, mode="w");
    h5createFile(file=poi.file);

    if(verbose) cat("added ind id into h5 file\n");
    #file.con[["individuals"]] <- ind.id
  # h5createGroup(file=poi.file, "individuals");
    h5write(ind.id, level=4, file=poi.file, name="individuals");

    if(verbose) cat("added poi id into h5 file\n");

    #file.con[["predictors_of_interest"]] <- poi.id
  #  h5createGroup(file=poi.file, "predictors_of_interest");
    h5write(poi.id, level=4, file=poi.file, name="predictors_of_interest");
  }
  num.poi.blocks <- ceiling(num.poi/poi.chunk.size);


  if(poi.type=="genotypes") {
    if(poi.file.type == "h5") {
      h5createDataset(
        file=poi.file, dataset="values", dims=c(num.ind,num.poi),
        storage.mode = "integer", chunk=c(num.ind,poi.chunk.size), 
        level=poi.compression.level
      )
    }

	for(j in 1:num.poi.blocks) {

      block.index <- ((j-1)*poi.chunk.size+1):(min(poi.chunk.size*j, num.poi));
      block.size <- length(block.index);
      maf <- runif(block.size, min=0.05, max=0.5);
      miss.rate <- runif(block.size, min=0, max=0.1);

      values <- matrix(0, ncol=block.size, nrow=num.ind);
      for(i in 1:block.size) {
        dosage.val <- runif(num.ind);
        geno.val <- integer(num.ind);
        geno.val[dosage.val < (1-maf[i])^2] <- 0;
        geno.val[((1-maf[i])^2 < dosage.val) & (dosage.val < (1 - maf[i]^2))] <- 1;
        geno.val[dosage.val > (1- maf[i]^2)] <- 2;
        geno.val[runif(num.ind)>1-miss.rate[i]] <- NA;
        values[,i] <- geno.val;
      }
      if(poi.file.type == "h5") {
        h5write(values, file=poi.file, name="values", index=list(NULL,block.index));
      }
      else {
        start <- ((j-1)*poi.chunk.size+1)
        end <- (min(poi.chunk.size*j, num.poi))
        values <- t(values)

        if (j == 1) {
          values <- cbind(poi.id[start:end], values)
          colnames(values) <- c("", ind.id)
          write.table(values, file=poi.txt.file, sep = "\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE)
        }
        else {
          values <- cbind(poi.id[start:end], values)
          write.table(values, file=poi.txt.file, sep = "\t", row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)
        }
      }

      if(verbose & (j %in% ceiling(seq(0.1,1,0.1)*num.poi.blocks))) cat("completed poi generation for ", j, " blocks out of ", num.poi.blocks, "\n");

	}
} else {
  if(poi.file.type == "h5") {
	  h5createDataset(file=poi.file, dataset="values", c(num.ind,num.poi),
                storage.mode = "double",
				chunk=c(num.ind,poi.chunk.size), level=poi.compression.level)
  }


	for(j in 1:num.poi.blocks) {
		block.index <- ((j-1)*poi.chunk.size+1):(min(poi.chunk.size*j, num.poi));
		block.size <- length(block.index);

		maf <- runif(block.size, min=0.05, max=0.5);
		miss.rate <- runif(block.size, min=0, max=0.1);

		values <- matrix(0.000, ncol=block.size, nrow=num.ind);
		for(i in 1:block.size) {
			geno.val <- ((1-maf[i])^2)*runif(num.ind,min=0,max=0.4)+2*(1-maf[i])*maf[i]*runif(num.ind,min=0.25, max=0.75)+(maf[i]^2)*runif(num.ind, min=0.6,max=1);
			geno.val[runif(num.ind)>1-miss.rate[i]] <- NA;
			values[,i] <- geno.val;
		}

    if(poi.file.type == "h5") {
        h5write(values, file=poi.file, name="values", index=list(NULL,block.index));
    } else {
      start <- ((j-1)*poi.chunk.size+1)
      end <- (min(poi.chunk.size*j, num.poi))
      values <- t(values)

      if (j == 1) {
        values <- cbind(poi.id[start:end], values)
        colnames(values) <- c("", ind.id)
        write.table(values, file=poi.txt.file, sep = "\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE)
      }
      else {
        values <- cbind(poi.id[start:end], values)
        write.table(values, file=poi.txt.file, sep = "\t", row.names=FALSE, col.names=FALSE, na="", quote=FALSE, append=TRUE)
      }
    }
		if(verbose & (j %in% ceiling(seq(0.1,1,0.1)*num.poi.blocks))) cat("completed poi generation for ", j, " blocks out of ", num.poi.blocks, "\n");
	}
}
  if(verbose) cat("added POI values into h5 file\n")

  #file.con[["values"]] <- t(values);
  #file.con$close_all();
  h5closeAll();

  invisible(list("num.poi"=num.poi, "num.ind"=num.ind, "num.resp.beta0"=num.resp.beta0, "num_resp.beta"=num_resp.beta, "num.resp.sd"=num.resp.sd, "bin.resp.beta0"=bin.resp.beta0, "bin_resp.beta"=bin_resp.beta));
}
