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
