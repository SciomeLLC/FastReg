
logisticRegression <- function(Y, G, X=NULL, Z=NULL, G.transform.fun=NULL, G.filter.fun=NULL, G.transform.args=NULL, G.filter.args=NULL, no.intercept=FALSE, max.iter=6, ...) {
  if(is.null(Z)) Z <- matrix(1, nrow=nrow(Y), ncol=1,dimnames=list(rownames(Y),"POI"));

  if(is.null(X)) {
      if(no.intercept) {
          X <- matrix(1, nrow=nrow(Y), ncol=0, dimnames=list(rownames(Y),NULL));
      } else {
          X <- matrix(1, nrow=nrow(Y), ncol=1, dimnames=list(rownames(Y), "Intercept"));
      }
  }

  common.subjects <- sort(intersect(rownames(Y), intersect(rownames(G), intersect(rownames(X), rownames(Z)))));

  X.index <- match(common.subjects, rownames(X));
  Y.index <- match(common.subjects, rownames(Y));
  G.index <- match(common.subjects, rownames(G));
  Z.index <- match(common.subjects, rownames(Z));

  G <- G[G.index, ,drop=FALSE];
  X <- X[X.index, ,drop=FALSE];
  Y <- Y[Y.index, ,drop=FALSE];
  Z <- Z[Z.index, ,drop=FALSE];

  remove(list=c("G.index", "X.index", "Y.index", "Z.index"));

  if(!is.null(G.filter.fun)) {
    filter.out <- do.call(G.filter.fun, c(list(G=G), G.filter.args));
    G <- G[,which(filter.out$keep==1),drop=FALSE]
  } else {
      filter.out <- NULL;
  }

  if(is.null(G.transform.args)) {
    G.transform.args <- list(effect.type="dosage");
  }

  if(!is.null(G.transform.fun)) {
    G <- do.call(G.transform.fun, c(list(G=G), G.transform.args));
  }


  xy.na <- which(rowSums(is.na(X))>0 |  rowSums(is.na(Z))>0 | is.na(Y));
  if(length(xy.na)>0) {
    X <- X[-xy.na, ,drop=FALSE];
    Y <- Y[-xy.na, ,drop=FALSE];
    G <- G[-xy.na, ,drop=FALSE];
    Z <- Z[-xy.na, ,drop=FALSE]
  }

  remove(list="xy.na");


  nS <- nrow(Y);
  nG <- ncol(G);
  nX <- ncol(X);
  nZ <- ncol(Z);

  POI_NAMES <- colnames(G);


  PARM_NAMES <- c(colnames(X), gsub("Intercept[:]","",paste0(colnames(Z), ":POI")));

  nP <- nZ+nX;


  se_beta <- beta <- matrix(0.0, ncol=nG, nrow=nP, dimnames=list(PARM_NAMES, POI_NAMES));

  XtX <- array(0.0, dim=c(nP, nP, nG), dimnames=list(PARM_NAMES, PARM_NAMES, POI_NAMES));
  XtY <- matrix(0.0, nrow=nP, ncol=nG, dimnames=list(PARM_NAMES, POI_NAMES));

  f_t <- matrix(1, ncol=nG, nrow=nS, dimnames=list(rownames(Y),POI_NAMES));
  for(v in 1:nG) {
      G.na <- which(is.na(G[,v]));
      f_t[G.na,v] <- 0;
      G[G.na,v] <- 0;
  }

  remove(list=c( "G.na"));

  One.V <- matrix(1, ncol=nG, nrow=1, dimnames=list(NULL, POI_NAMES));

  ################################
  start.time <- proc.time()

  max.iter <- max.iter + 1;

  for(iter in 1:max.iter) {

	  mu <- X%*%beta[1:nX, , drop=FALSE]; ## NUM_SUBJECTS x NUM_VARIANT
    for(j in 1:nZ) mu <- mu + G*(Z[,j,drop=FALSE]%*%beta[nX + j,,drop=FALSE]);
    mu <- 1/(1+exp(-mu))
	  W_t <-  f_t * mu * (1-mu)

    ## Compute Residuals z and
    ## then compute X'z
    res  <- (Y%*%One.V-mu);

    ### Set residuals and weight observations with missing variant call to 0;

	  remove(list="mu");

    ################ Using custom written function because 3D array multiplication is not available in R


	  XtX[1:nX,1:nX,] <- matmult3D(W_t, X);

    for(j in 1:nZ) {
      XtX[1:nX,nX+j,] <- matmult2D(W_t, X, G*Z[,j]);
      for(d in 1:nX) XtX[nX+j, d,] <- XtX[d,nX+j,];
	    for(k in 1:nZ) XtX[nX+j,nX+k,]   <- matmult1D(W_t, G*Z[,j],G*Z[,k]);
   }


	if(iter!=max.iter) {
		XtY[1:nX,] <- matmult2D(f_t, X, res);
		for(j in 1:nZ) XtY[nX+j,] <- matmult1D(f_t, res, G*Z[,j]);
    beta <- beta  + solve3D(a=XtX, b=XtY);
	} else {
		se_beta <- sqrt(abs(diag3D(solve3D(a=XtX))));
	}
}

N <- colSums(f_t);
DF <- N-nP;

end.time <- proc.time();

elapsed.time <- (end.time-start.time);

#cat("Successful Execution Took with Elapsed Time:", elapsed.time, "\n");

return(list("beta"=beta, "se_beta"=se_beta, "N"=N,"DF"=DF, "filter.out"=filter.out, "elapsed.time"=elapsed.time));

}



