rrr_big <- function(X, Y=NULL, power=1, K.X=min(dim(X)) - 1, K.Y=1,verbose=T,lr.file = tempfile(),ncores=1,lr.return=F) {

  if (is.null(Y)) {
    if (all(vapply(c('matrix','array','FBM'), function(x) inherits(X,x),TRUE) == FALSE)) 
      stop("X has to be an object of class 'array', 'matrix' or 'FBM'.")
    if (K.X > min(dim(X))| floor(K.X) != K.X) stop("K.X has to be and integer lower than or equal to min(dim(X)).")
  }
  else {
    if (all(vapply(c('matrix','array','FBM'), function(x) inherits(X,x),TRUE) == FALSE)) 
      stop("X has to be an object of class 'matrix' or 'FBM'.")
    
    if (all(vapply(c('matrix','array','FBM'), function(x) inherits(Y,x),TRUE) == FALSE)) stop("Y has to be an object of class 'matrix' or 'FBM'.")
    stop("Y have to be an object of class 'array', 'matrix' or 'FBM'.")
    
    if (nrow(Y) != nrow(X)) stop("Both matrices need to have the same number of rows.")
    if (K.X > min(dim(X))| floor(K.X) != K.X) stop("K.X has to be an integer lower than or equal to min(dim(X)).")
    if (K.Y > min(dim(X),dim(Y)) | floor(K.Y) != K.Y) stop("K.Y has to be an integers lower than or equal to min([dim(X),dim(Y)]).")
  }
  if (!any(power %in% c(1,-1))) stop("power has to be '1' or '-1'")

  #Getting dimensions.
  n <- nrow(X); p_Y <- ncol(Y); p_X <- ncol(X)

  #Getting SVD of X.
  if (verbose) message("Getting SVD of predictors matrix (dimension ",n," x ",p_X , ").")
  if (inherits(X, "FBM") == TRUE) SVD.x <- bigstatsr::big_randomSVD(X, k = K.X, verbose = verbose,ncores = ncores)
  else {
    SVD.x <- svd(X, nu = K.X, nv=K.X)
    SVD.x$d <- SVD.x$d[1:K.X]
  }

  #Retaining axes of positive eigen-values to avoid singularity.
  tmp <- SVD.x$d > 1e-5
  if (sum(tmp) < min(n, p_X) - 1) {
    SVD.x$u <- SVD.x$u[, tmp]
    SVD.x$v <- SVD.x$v[, tmp]
    SVD.x$d <- SVD.x$d[tmp]
  }
  if (length(tmp) == 0) stop("No positive eigenvalues for X. Maybe a larger K?")
  if (is.null(Y)) return(SVD=SVD.x)
  else {
    #Creating matrix R.
    if (verbose) message("Creating matrix 'R' (dimension ",K.X," x ",p_Y , ").")
    if (inherits(Y, "FBM")) {
      R <- bigstatsr::FBM(K.X, p_Y, create_bk = T)$save()
      #Filling up matrix R.
      count <- 0
      bigstatsr::big_apply(Y, a.FUN = function(y, ind) {
        if (verbose) cat("Working chunk ",(count <<- count + 1),"of",ceiling(p_Y / bigstatsr::block_size(p_Y, 1)),".\n")
        R[, ind] <- crossprod(SVD.x$u, y[,ind] * p_Y)
        NULL
      },
      a.combine = 'c',
      ind = seq_len(p_Y),
      block.size = bigstatsr::block_size(p_Y, 1))
    }
    else R <-  crossprod(SVD.x$u, Y * p_Y)

    #Creating matrix L.
    if (verbose) message("Creating matrix 'L' (dimension ",p_X," x ",K.X , ").")
    if (inherits(X, "FBM") == TRUE) L <- bigstatsr::as_FBM(SVD.x$v /  p_X)
    else L <- SVD.x$v /  p_X
    for (k in 1 : K.X) L[, k] <- L[, k] * (SVD.x$d[k] ^ power)

    #Creating matrix B.
    if (verbose) message("Getting 'L*R' (dimensions ",p_X," x ",p_Y , ").")
    if (inherits(L, "FBM") == TRUE | inherits(R, "FBM") == TRUE) {
      LR <- bigstatsr::FBM(p_X, p_Y, create_bk = T,backingfile = lr.file)$save()
      count <- 0
      if (inherits(L, "FBM") == TRUE) {
        bigstatsr::big_apply(LR, a.FUN = function(y, ind) {
          if (verbose) cat("Working chunk ",(count <<- count + 1),"of",ceiling(p_Y / bigstatsr::block_size(p_Y, 1)),".\n")
          LR[, ind] <- bigstatsr::big_prodMat(L, R[,ind], block.size = bigstatsr::block_size(K.X, 1))
          NULL
        },
        a.combine = 'c',
        ind = seq_len(p_Y),
        block.size = bigstatsr::block_size(p_Y, 1))
      }
      else {
        bigstatsr::big_apply(LR, a.FUN = function(y, ind) {
          if (verbose) cat("Working chunk ",(count <<- count + 1),"of",ceiling(p_Y / bigstatsr::block_size(p_Y, 1)),".\n")
          LR[, ind] <- L %*% R[,ind]
          NULL
        },
        a.combine = 'c',
        ind = seq_len(p_Y),
        block.size = bigstatsr::block_size(p_Y, 1))
      }
    }
    else LR <- L %*% R
    if (lr.return) return(LR)
    else {
      if (verbose) message("Getting SVD of LR.")
      if (inherits(LR, "FBM") == TRUE) SVD <- bigstatsr::big_randomSVD(LR, k = K.Y, verbose = verbose,ncores=ncores)
      else SVD <- svd(LR, nu = K.Y, nv=K.Y)
      return(list("SVD"=SVD,"LR"=LR))
    }
  }
}

