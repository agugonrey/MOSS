prepro_sub <- function(X,scale.,norm.) {
  if(any(vapply(c("matrix","FBM","array"), function (x) inherits(X,x),TRUE))) {
    if (inherits(X, "FBM") == TRUE) {
        bigstatsr::big_apply(X, function(x, ind) {
        if (scale.) X[,ind] <- scale(x[,ind])
        if (norm.) X[,ind] <- x[,ind] / sqrt(ncol(X))
      },
      a.combine = 'c',
      ind = seq_len(ncol(X)),
      block.size = bigstatsr::block_size(ncol(X), 1))
    }
    else {
      if (scale.) X <- scale(X)
      if (norm.) X <- X / ncol(X)
    }
  }
  else stop("Only objects of class 'array', 'matrix' or 'FBM' supported.")
  return(X)
}

