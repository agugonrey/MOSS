#' Missing values imputation by the mean of each column.
#'
#' This function is called by moss to count the impute missing values
#' by the mean of each column within omic blocks.
#' If any column has more than 20% missing values an error is displayed.
#'
#' Meant for objects of class 'matrix', 'FBM', or 'array'. 
#' @param X An object of class 'matrix', 'FBM', or 'array'.
#' @return Returns input with imputed missing values.
#' @export
prepro_na <- function(X) {
  if (any(vapply(
    c("matrix", "FBM", "array"),
    function(x) inherits(X, x), TRUE
  ))) {
    n <- nrow(X)
    if (inherits(X, "FBM") == TRUE) {
      prop_na <- numeric(ncol(X))
      bigstatsr::big_apply(X, function(X, ind) {
        # access a subset of columns as a standard R matrix
        X.sub <- X[, ind, drop = FALSE]
        # get the location (i, j) of missing values
        ind_na <- which(is.na(X.sub), arr.ind = TRUE)
        # Check that there are not more than 20% NA.
        if(any(colSums(is.na(X.sub[,unique(ind_na[,2])])) / n > 0.2)) 
          stop("At least one features with more than 20% missing!")
        # compute the corresponding mean for each column j
        means <- colMeans(X.sub, na.rm = TRUE)[ind_na[, 2]]
        # update j (relative to subset) to global 'ind'
        ind_na[, 2] <- ind[ind_na[, 2]]
        # fill positions with corresponding means
        X[ind_na] <- means
        # here we don't want to return anything, so `NULL`.
        NULL
      }, a.combine = 'c')
    }
    else {
      X <- apply(X, 2, function(x) {
        tmp <- is.na(x)
        if(sum(tmp) / n > 0.2) 
          stop("At least one features with more than 20% missing!")
        x[tmp] <- mean(x,na.rm = T)
        x
      })
    }
  }
  else {
    stop("Only objects of class 'array', 'matrix' or 'FBM' supported.")
  }
  return(X)
}
