#' Number of missing values within a matrix.
#' 
#' This function is called by moss to count the number of missing values within the (extended) matrices.
#' 
#' Ment for objects of class 'matrix', 'FBM', or 'array'.
#' @param X An object of class 'matrix', 'FBM', or 'array'.
#' @return Returns total number of missing values in X.
#' @export
prepro_na <- function(X) {
  if(any(vapply(c("matrix","FBM","array"), function (x) inherits(X,x),TRUE))) {
    if (inherits(X, "FBM") == TRUE) {
      na_count <- sum(bigstatsr::big_apply(X, function(x, ind) {
       sum(is.na(x[,ind]))
      },
      a.combine = 'c',
      ind = seq_len(ncol(X)),
      block.size = bigstatsr::block_size(ncol(X), 1)))
    }
    else {
      na_count <- sum(apply(X,2,function(x) sum(is.na(x))))
    }
  }
  else stop("Only objects of class 'array', 'matrix' or 'FBM' supported.")
  return(na_count)
}
