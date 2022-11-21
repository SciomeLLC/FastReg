#'This file consist of trivial 3D array operation functions not intended to be used directly by end-user

# matrix multiplication of 3D matrices
matmult3D <- function(w, a, b=NULL) {
  p <- dim(w)[2];
  m <- dim(a)[2];


  if(is.null(b)) {
    v <- array(0, dim=c(m, m, p), dimnames=list(dimnames(a)[[2]], dimnames(a)[[2]], dimnames(w)[[2]]));
    for(d in 1:p) v[,,d] <- t(a)%*%(a*w[,d]);

  }  else {
    n <- dim(b)[2];
    v <- array(0, dim=c(m, n,p), dimnames=list(dimnames(a)[[2]], dimnames(b)[[2]], dimnames(w)[[2]]));
    for(d in 1:p) v[,,d] <- t(a)%*%(b*w[,d]);
  }

  return(v)
}

# matrix multiplication of 2D matrices
matmult2D <- function(w, a, b=NULL) {

  if(is.null(b)) {
      v <- t(a)%*%(w*a);
  } else {
    v <- t(a)%*%(w*b);
  }

  return(v)
}

# matrix multiplication of 1D matrices
matmult1D <- function(w, a, b=NULL) {

  if(is.null(b)){
    v <- colSums(w*a*a);
  } else {
    v <- colSums(w*a*b)
  }

  return(v)
}

##Generalized inverse function imported from MASS package
ginv <- function (X, tol = sqrt(.Machine$double.eps))  {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}


gsolve3D <- function(a,b=NULL) {
  p <- dim(a)[3];
  if(!is.null(b)) {
    i <- matrix(0, nrow=dim(a)[1], ncol=dim(a)[3],dimnames=list(dimnames(a)[[1]], dimnames(a)[[3]]));
    for(d in 1:p) {
      i[,d] <- ginv(a[,,d])%*%b[,d];
    }
  } else {
    i <- array(0, dim=dim(a), dimnames=dimnames(a));
    for(d in 1:p) {
      i[,,d] <- ginv(a[,,d]);
    }
  }
  return(i)
}


solve3D <- function(a, b=NULL) {
  p <- dim(a)[3];
  if(!is.null(b)) {
    i <- matrix(0, nrow=dim(a)[1], ncol=dim(a)[3],dimnames=list(dimnames(a)[[1]], dimnames(a)[[3]]));
    for(d in 1:p) {
      i[,d] <- solve(a=a[,,d],b=b[,d]);
    }
  } else {
    i <- array(0, dim=dim(a), dimnames=dimnames(a));
    for(d in 1:p) {
      i[,,d] <- solve(a=a[,,d]);
    }
  }
  return(i)
}

matmultSolve <- function(a,b) {
  p <- dim(a)[3];
  i <- matrix(0, nrow=dim(a)[1], ncol=dim(a)[3], dimnames=list(dimnames(a)[[1]], dimnames(a)[[3]]));
  for( d in 1:p) {
    i[,d] <- a[,,d]%*%b[,d]
  }
  return(i)
}

matmultElement <- function(a,b) {
  p <- dim(a)[3];
  for(d in 1:p) {
    a[,,d] <- a[,,d]*b[d]
  }
  return(a)
}

# returns the diagonal of a matrix
diag3D <- function(a) {
  p <- dim(a)[3];
  v <- matrix(0, nrow=dim(a)[1], ncol=dim(a)[3], dimnames=list(dimnames(a)[[1]], dimnames(a)[[3]]));
  for(d in 1:p) {
    v[,d] <- diag(a[,,d]);
  }
  return(v);
}



