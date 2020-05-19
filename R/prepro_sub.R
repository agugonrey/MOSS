#' Scale and normalize columns of a matrix.
#' 
#' This function is called by moss to scale and normalize (extended) matrices.
#' 
#' Ment for objects of class 'matrix', 'FBM', or 'array'.
#' @param X An object of class 'matrix', 'FBM', or 'array'.
#' @param scale.arg Should we scale columns? Logical.
#' @param norm.arg Should we normalize columns? Logical.
#' @return A matrix with scaled and/or normalized columns.
#' @export
prepro_sub <- function(X,scale.arg,norm.arg) {
  if(any(vapply(c("matrix","FBM","array"), function (x) inherits(X,x),TRUE))) {
    if (inherits(X, "FBM") == TRUE) {
        bigstatsr::big_apply(X, function(x, ind) {
        if (scale.arg) X[,ind] <- scale(x[,ind])
        if (norm.arg) X[,ind] <- x[,ind] / sqrt(ncol(X))
      },
      a.combine = 'c',
      ind = seq_len(ncol(X)),
      block.size = bigstatsr::block_size(ncol(X), 1))
    }
    else {
      if (scale.arg) X <- scale(X)
      if (norm.arg) X <- X / ncol(X)
    }
  }
  else stop("Only objects of class 'array', 'matrix' or 'FBM' supported.")
  return(X)
}

