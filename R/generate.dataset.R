#' generate.dataset function to generate Simulated Dataset. NOTE THAT num.subj = num.subj saved + 1
#'@param seed random number generator seed (default=22342)
#'@param num.subj number of subjects, default = 10000;
#'@param num.variants number of variants, default = 500;
#'@param num.pc number of covariates excludinf Intercept, Age, Sex, default=5;
#'@param num.pheno number of phenotypes to simulate, default=2;
#'@param add.intercept logical denoting whether add intercept column in X matrix
#'@param output.dir output directory to save generate RDS objects;
#'@param prefix character string used as prefix for RDS filenames;
#'@return 1 if successful
#'@export
#'@importFrom data.table fwrite
#'@author: Deepak Mav (deepak.mav@sciome.com)
#'Last Edit Date: 07/26/2024

generate.dataset <- function(seed=22342, num.pc=5, num.subj=10000, num.variants=500, num.pheno=2, add.intercept = FALSE, output.dir=".", prefix="") {
  set.seed(seed);
  
  ### Generate Covariate Matrix (X);
  PC.names <- paste0("PC",formatC(1:num.pc, format="d", flag="0", digits=floor(log10(num.pc))));
  
  
  X <- matrix(NA, nrow=num.subj,ncol=num.pc+4, dimnames=list(NULL, c("Age", "Sex", "AgeSq", "AgeSqTimesSex", PC.names)));
  X[, "Age"] <- pmin(pmax(1,rnorm(num.subj, mean=35, sd=5)), 85);
  X[,"AgeSq"] <- X[, "Age"]^2;
  X[,"Sex"] <- 1*(runif(num.subj)>0.5);
  X[, "AgeSqTimesSex"] <- X[,"AgeSq"]*X[,"Sex"];
  
  for(j in 1:num.pc) X[, PC.names[j]] <- rnorm(num.subj, mean=rnorm(1,0,1), sd=runif(1, 0.5, 1.5));
  
  
  ### Generate Genotype Call Matrix (X); 
  maf <- runif(num.variants, 0.001, 1);		   	
  G <- do.call(cbind, lapply(maf, function(x) rbinom(num.subj, size=2, prob=x)));
  colnames(G) <- paste0("v", formatC(1:num.variants, format="d", flag="0", digits=floor(log10(num.variants))));
  
  pheno.names <- paste0("Y",formatC(1:num.pheno, format="d", flag="0", digits=floor(log10(num.pheno))));
  
  Y <- matrix(NA, nrow=num.subj, ncol=num.pheno, dimnames=list(NULL, pheno.names));
  
  q <- ncol(X);  
  
  ### Generate Phenotype Matrix;
  for(i in 1:num.pheno) {
    beta.X <- matrix(runif(q,min=-0.25,max=0.25), nrow=q, ncol=1);
    beta.G <- matrix(rnorm(num.variants, mean=0.1, sd=0.05), nrow=num.variants, ncol=1);
    Y[,i] <- rnorm(num.subj, mean=rpois(1,20)+X%*%beta.X+G%*%beta.G, sd=runif(1));
  }
  
  
  ### Generate (1%) random missingness in Y, X and G; 
  for(i in 1:num.pheno) Y[runif(num.subj)<0.01,i] <- NA;
  for(i in 1:q) X[runif(num.subj)<0.01,i] <- NA;
  for(i in 1:num.variants) G[runif(num.subj)<0.01,i] <- NA;
  
  if(add.intercept) X <- cbind("Intercept"=rep(1, nrow=num.subj), X);
  
  # G = rbind(G,rbinom(num.variants, size=2, prob=maf[1]));
  X = X[-nrow(X),]
  Y = Y[-nrow(Y),]
  dir.create(output.dir, showWarnings=FALSE);
  saveRDS(Y, file=file.path(output.dir, paste0(prefix, "Y.rds")));
  saveRDS(X, file=file.path(output.dir, paste0(prefix, "X.rds")));
  # saveRDS(G, file=file.path(output.dir, paste0(prefix, "G.rds")));
  rm(X,Y)
  
  poi.id = colnames(G);
  ind.id = paste0("IND",1:num.subj)
  
  num.poi <- length(poi.id)
  num.ind <- length(ind.id)
  
  seqs <- c(
    "chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559, "chr4" = 190214555, "chr5" = 181538259,
    "chr6" = 170805979, "chr7" = 159345973, "chr8" = 145138636, "chr9" = 138394717, "chr10" = 133797422,
    "chr11" = 135086622, "chr12" = 133275309, "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189,
    "chr16" = 90338345, "chr17" = 83257441, "chr18" = 80373285, "chr19" = 58617616, "chr20" = 64444167,
    "chr21" = 46709983, "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415, "chrMT" = 16569
  )
  
  seq.names <- names(seqs)
  seq.names <- gsub("chr", "", seq.names)
  num.seq <- length(seqs)
  
  cat("Generating .bim file\n")
  poi.seq <- as.integer(cut(sort(runif(num.poi)), breaks = c(0, cumsum(as.numeric(seqs / sum(as.numeric(seqs)))))))
  poi.pos <- integer(num.poi)
  poi.ref <- poi.alt <- character(num.poi)
  for (seq in 1:num.seq) {
    w <- which(poi.seq == seq)
    num.w <- length(w)
    if (num.w > 0) {
      poi.pos[w] <- sort((unique(ceiling(seqs[seq] * runif(num.w * 2))))[1:num.w])
      rf <- ceiling(4 * runif(num.w))
      poi.ref[w] <- c("A", "T", "C", "G")[rf]
      poi.alt[w] <- c("A", "T", "C", "G")[4 - rf + 1]
    }
  }
  
  poi.df <- as.data.table(data.frame(CHROM = seq.names[poi.seq], ID = poi.id, POS = integer(num.poi), "COORD" = poi.pos, "ALT" = poi.alt, "REF" = poi.ref, stringsAsFactors = FALSE, check.names = FALSE))
  remove(list = c("poi.seq", "poi.pos", "poi.ref", "poi.alt", "w", "num.w", "rf", "seq", "seqs", "num.seq", "seq.names"))
  write.table(poi.df, file = file.path(output.dir, paste0(prefix,"p.bim")), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "")
  remove(list = "poi.df")
  
  cat("Generating .fam file\n")
  
  fam <- data.frame(
    `#FID` = rev(ind.id), 
    IID = ind.id,
    check.names = FALSE
  )
  
  write.table(fam, file = file.path(output.dir, paste0(prefix,"p.fam")), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE, na = "")
  remove(list = "fam")
  
  
  cat("Generating .bed file\n")
  outfile = file.path(output.dir, paste0(prefix,"p.bed"))
  if (file.exists(outfile)) unlink(outfile)
  con <- file(outfile, "wb")
  poi.chunks <- 1
  nBytes <- ceiling(num.ind / 4)
  
  d <- matrix(c(
    0, 1, 1, 0, 1, 1, 0, 0,
    0, 0, 0, 1, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1
  ), nrow = 8)
  
  xd <- as.raw(drop(matrix(2^(7:0), nrow = 1, ncol = 8) %*% d))
  temp_output <- writeBin(xd, con = con)
  
  base.vec <- matrix(2^(0:7), nrow = 1, ncol = 8)
  
  poi.index <- 1:num.poi
  chunk.size <- length(poi.index)
  poi.chunk <- G
  
  xd <- matrix(0, nrow = nBytes, ncol = chunk.size)
  
  for (k in 1:chunk.size) {
    d <- matrix(0, nrow = 8, ncol = nBytes)
    ip <- which(poi.chunk[, k] == 1)
    if (length(ip) > 0) {
      d[2 * (ip - 1) + 1] <- 0
      d[2 * (ip - 1) + 2] <- 1
    }
    ip <- which(poi.chunk[, k] == 0)
    if (length(ip) > 0) {
      d[2 * (ip - 1) + 1] <- 1
      d[2 * (ip - 1) + 2] <- 1
    }
    ip <- which(is.na(poi.chunk[, k]))
    if (length(ip) > 0) {
      d[2 * (ip - 1) + 1] <- 1
      d[2 * (ip - 1) + 2] <- 0
    }
    xd[, k] <- drop(base.vec %*% d)
    remove(list = "d")
  }
  
  # print(xd[, k])
  data <- as.raw(c(xd))
  # print(data)
  writeBin(data, con=con); flush(con);
    
  remove(list = c("poi.chunk", "xd"))
  
  close(con)
  
  invisible(1)
}
