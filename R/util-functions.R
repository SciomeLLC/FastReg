POI.filter <- function(G, maf.threshold=0.01, hwe.threshold=0.05) {


  nS <- nrow(G)-colSums(1*is.na(G));
  a.freq <- colSums(G, na.rm=TRUE)/(2*nS);
  b.freq <- colSums(2-G, na.rm=TRUE)/(2*nS);

  maf.freq <- pmin(a.freq, b.freq);
  aa.of <- colSums(G==0, na.rm=TRUE);
  ab.of <- colSums(G==1, na.rm=TRUE);
  bb.of <- colSums(G==2, na.rm=TRUE);
  p <- (2*aa.of + 1*ab.of)/(2*nS);
  aa.ef <- (p^2)*nS
  ab.ef <- 2*p*(1-p)*nS;
  bb.ef <- (1-p)^2*nS;

  HWE.chisq <- (aa.of-aa.ef)^2/aa.ef + (ab.of-ab.ef)^2/ab.ef + (bb.of-bb.ef)^2/bb.ef;
  HWE.pval <- 1-pchisq(HWE.chisq, df=1);
  keep <- 1*((maf.freq>=maf.threshold) & (HWE.pval>=hwe.threshold))

  data.frame("POI"=colnames(G), "Freq (a)"=a.freq, "Freq (b)"=b.freq, "MAF"=maf.freq, "HWE Chisq"=HWE.chisq, "HWE Pvalue"=HWE.pval, "keep"=keep,
                check.names=FALSE, stringsAsFactors=FALSE);

}

POI.transform <- function(G, effect.type="additive") {
  if(!(effect.type %in% c("dosage", "additive", "dominant", "recessive"))) stop("invalid effect.type specified");


  if(effect.type == "dominant") {
      G[G==2] <- 1;
  } else if(effect.type == "recessive") {
      G[G==1] <- 0;
      G[G==2] <- 1;
  }

  return(G)
}


add.covar <- function(df, cv, cv.type, standardize=FALSE, X.mat=NULL, levels=NULL, ref.level=NULL, colinearity.rsq=0.99, verbose=1) {
  nS <- nrow(df);
  if(is.null(X.mat)) X.mat <- matrix(0.0, nrow=nS, ncol=0, dimnames=list(rownames(df), NULL));
  if(nrow(X.mat)!=nS) stop("number of rows for X.mat and df argument do not match")

  if(!(cv %in% colnames(df))) {
    stop(paste(cv,  "is not present in df"));
  }

  if(!cv.type %in% c("numeric", "categorical")) {
    stop("cv.type must be either numeric, categorical");
  }

  if((cv.type == "numeric") & inherits(df[,cv], "character")) {
    stop("numeric cv.type for non-numeric variable used");
  }

  candidate.cols <- list();
  if(cv.type=="numeric") {
    candidate.cols[[cv]] <- df[,cv,drop=TRUE];
  } else {
      if(is.null(levels)) levels <- sort(unique(df[,cv,drop=TRUE]));
      if(is.null(ref.level)) ref.level <- levels[1];
      if(!(ref.level %in% levels)) stop("ref.level must be one of unique value of var");

      levels <- setdiff(levels,ref.level) ;
      if(ncol(X.mat)==0) levels <- c(ref.level,levels);

      nl <- length(levels);
      for(lev in levels) candidate.cols[[lev]] <- 1*(df[,cv] == lev);

      if(ncol(X.mat)==0) {
        names(candidate.cols) <- paste0(cv, "(", names(candidate.cols),")");
      } else {
        names(candidate.cols) <- paste0(cv, "(", names(candidate.cols), " vs. ", ref.level,")");
      }
      if(standardize) stop("standardization of non-numeric variable not permitted")
  }

  if(standardize) {
    candidate.cols <- lapply(candidate.cols, function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=FALSE));
  }

  nl <- length(candidate.cols);
  nS <- nrow(df);
  retain <- rep(FALSE,nl);



  for(d in 1:nl) {

    if(ncol(X.mat)==0 & d==1) {
        retain[d] <- TRUE;
        next;
    }

    y <- matrix(candidate.cols[[d]], nrow=nS, ncol=1, dimnames=list(rownames(df),NULL));
    Z <- do.call(cbind, c(list(X.mat), candidate.cols[retain]));
    w <- which(!is.na(y) & rowSums(is.na(Z))==0);
    Z <- Z[w,];
    y <- y[w];

    ssa <- sum(y^2);
    sse<- sum((y-Z%*%ginv(t(Z)%*%Z)%*%(t(Z)%*%y))^2);
    rsquared <- 1-sse/ssa;
    remove(list=c("y","Z","w"))
    if(rsquared<=colinearity.rsq) {
        retain[d] <- TRUE;
    } else {
      cat(candidate.cols[d], "was not added to design matrix due to potential of colinearity\n");
    }
  }

  candidate.cols <- candidate.cols[which(retain)];
  Z <- do.call(cbind, c(list(X.mat), candidate.cols[retain]));

  return(Z)

}


createDesignMatrix <- function(df, covariates=NULL, covariate.type=NULL, covariate.standardize=NULL, no.intercept=FALSE, covariate.levels=NULL, covariate.ref.level=NULL, colinearity.rsq=0.99, verbose=1) {
  if(no.intercept!=1) {
      X.mat <- matrix(1, nrow=nrow(df), ncol=1, dimnames=list(rownames(df), "Intercept"));
  } else {
      X.mat <- matrix(0, nrow=nrow(df), ncol=0, dimnames=list(rownames(df), NULL));
  }

  if(is.null(covariates)) return(X.mat);
  if(!all(covariates %in% colnames(df))) stop("invalid covariates argument")

  if(is.null(covariate.type)) {
    covariate.type <- sapply(covariates, function(cv) ifelse(inherits(df[,cv,drop=TRUE], "factor") | inherits(df[,cv,drop=TRUE], "character"), "numeric", "categorical"));
    names(covariate.type) <- covariates;
  }

  if(!identical(names(covariate.type), covariates)) stop("covariate.type names do not match covariates");
  if(!all(covariate.type %in% c("numeric", "categorical",  "count"))) stop("invalid covariate.type");

  ccov <- covariates[covariate.type %in% c("categorical")];


  if(is.null(covariate.levels)) {
     covariate.levels <- lapply(covariates, function(x) c());
     names(covariate.levels) <- covariates;
     for(cv in ccov) covariate.levels[[cv]] <- sort(unique(df[,cv,drop=TRUE]));
  }


  if(is.null(covariate.ref.level)) {
    covariate.ref.level <- lapply(covariates, function(x) c());
    names(covariate.ref.level) <- covariates;
    for(cv in ccov) covariate.ref.level[[cv]] <- covariate.levels[[cv]][1];
  }

  if(is.null(covariate.standardize)) {
      covariate.standardize <- rep(FALSE,length(covariates));
      names(covariates.standardize) <- covariates;
  }

  if(inherits(covariate.standardize, "character")) {
    covariate.standarize <- (covariates %in% covariate.standardize)
    names(covariates.standardize) <- covariates;
  }

  if(length(ccov) > 0) {
    for(cv in ccov) {
      wc <- which(!is.na(df[,cv,drop=TRUE]));
      if(!all(df[wc,cv,drop=TRUE] %in% covariate.levels[[cv]])) stop("invalid covariate.levels")
      if(!(covariate.ref.level[[cv]] %in% covariate.levels[[cv]])) stop("invalid covariate.ref.level")
    }
  }

  for(cv in covariates) {
    if(verbose>0) cat("adding", cv, "to design matrix\n")
    X.mat <- add.covar(df=df, cv=cv, cv.type=covariate.type[cv], standardize=covariate.standardize[cv], X.mat=X.mat, levels=covariate.levels[[cv]], ref.level=covariate.ref.level[[cv]], colinearity.rsq=colinearity.rsq);
  }

  return(X.mat);
}

convertPOI.TxtToH5 <- function(file, outfile, ...) {
  invisible(1)
}

open.POI.file <-function(type, file, n=NULL, p=NULL) {
  if(type=="h5") {
    object <- h5file(file, mode="r");
    # object <- h5file(file, mode="r");
  } else {
    object <- BEDMatrix(file, n=n, p=p)
  }
  return(object);
}

close.POI.file <- function(type, object) {
  if(type=="h5") {
    object$close_all();
  }
  invisible(1)
}

get.POI.names <- function(type, file) {

  if(type=="h5"){
      object <- h5file(file, mode="r");
      POI.names <- object[["predictors_of_interest"]][];
      object$close_all();
  } else  {

      bim.file <- gsub("[.]bed", ".bim", file)
      if(!file.exists(bim.file)) stop("bim file missing")
      bim.dt <- fread(file=bim.file, header=FALSE, sep="\t");
      colnames(bim.dt) <- c("chr", "vid", "mpos", "pos", "A", "B");

      POI.names <- paste0(bim.dt[,"chr",drop=TRUE], "_", bim.dt[,"pos",drop=TRUE], "_",bim.dt[,"A",drop=TRUE],"/", bim.dt[,"B",drop=TRUE])
  }
  return(POI.names);
}

get.POI.individuals <- function(type, file) {
  if(type=="h5") {
        object <- h5file(file, mode="r")
        POI.individuals <- object[["individuals"]][];
        object$close_all();
  } else {
      fam.file <- gsub("[.]bed", ".fam", file)
      if(!file.exists(fam.file)) stop("bim file missing")
      fam.dt <- fread(file=bim.file, header=FALSE, sep="\t");
      colnames(fam.dt)[1:5] <- c("FID", "IID", "PID", "MID", "sex");
      POI.individuals <- paste0(fam.dt[,"FID",drop=TRUE], "_", fam.dt[,"IID",drop=TRUE])

  }
  return(POI.individuals)
}


get.POI.Matrix <- function(type, file.object, POI.names=NULL, POI.individuals=NULL) {
  if(type=="h5") {

    if(is.null(POI.names)) POI.names <- file.object[["predictors_of_interest"]];
    POI.order <- order(POI.names);
    POI.index <- (1:length(POI.names))[POI.order];
    POI.names <- POI.names[POI.order];

    if(is.null(POI.individuals)) POI.individuals  <- file.object[["individuals"]];
    IND.order <- order(POI.individuals);
    IND.index <- (1:length(POI.individuals))[IND.order];
    IND <- POI.individuals[IND.order];

    function(ind, poi, file=file.object, POI=POI.names, POI.rd=POI.index, IND=POI.individuals, IND.rd=IND.index) {
      if(!all(ind %in% IND)) stop("invalid individual index");
      if(!all(poi %in% POI)) stop("invalid predictor_of_interest index");
      G <- file[["values"]][POI.rd[match(poi, POI)],IND.rd[match(ind, IND)]];
      dimnames(G) <- list(poi, ind);
      return(t(G));
    }
  } else {
    function(ind, poi, file=file.object, POI=POI.names, POI.rd=POI.index, IND=POI.individuals, IND.rd=IND.index) {
      if(!all(ind %in% IND)) stop("invalid individual index");
      if(!all(poi %in% POI)) stop("invalid predictor_of_interest index");
      G <- file[Ind.rd[match(ind, IND)], POI.rd[match(poi, POI)],drop=FALSE];
      dimnames(G) <- list(ind,poi);
      return(G);
    }
  }
}



