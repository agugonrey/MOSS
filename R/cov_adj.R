#' Adjust omic blocks for covariates effects.
#'
#' This function is called by moss to adjust a series of omic 
#' blocks for covariates effects.
#' @param data.blocks List containing omic blocks of class 'matrix' or
#' 'FBM'. In each block, rows represent subjects and columns features.
#' @param covs Covariates which effect we wish to adjust for. 
#' Accepts objects of class matrix, data.frame, numeric, or 
#' character vectors.
#' @param n Number of subjects. Numeric.
#' @param dim.names list of vectors with samples names, and features names
#' by omic block. If NULL, a list of artificial names is created.
#' Defaults to NULL.
#' @return Returns the covariates-adjusted elements in data.blocks.
#' @export
#' @examples
#' library("MOSS")
#' sim_data <- simulate_data()
#' set.seed(43)
#'
#' # Extracting simulated omic blocks.
#' sim_blocks <- sim_data$sim_blocks[-4]
#' 
#' # Using fourth block as covariates.
#' covs <- sim_data$sim_blocks[[4]]
#' 
#' # Adjust omic blocks for covariates effects.
#' sim_blocks_adj <- cov_adj(sim_blocks,covs,nrow(covs))

cov_adj <- function (data.blocks, covs, n, dim.names = NULL)
{
  M <- length(data.blocks)
  crossprod.na <- function(x,y=x) {
    crossprod(replace(x,is.na(x),0),
              replace(y,is.na(y),0))
  }
  if (Reduce("+", lapply(c("matrix", "character", "numeric",
                           "array", "data.frame"), inherits, x = covs)) == 0) {
    stop("'covs' must be a vector, matrix, array, or dataframe")
  }
  if (is.null(dim.names)) {
    tmp_row_names <- seq_len(nrow(data.blocks[[1]]))
    dim.names <- lapply(data.blocks, function(x) list(tmp_row_names,
                                                      seq_len(ncol(x))))
  }
  if (is.null(dim(covs)))
    covs <- as.matrix(covs)
  if (nrow(covs) != n)
    stop("Different number or rows between covariates and omics")
  if (dim(covs)[2] > 1 & is.null(colnames(covs))) {
    warning("Missing column names in covs")
    colnames(covs) <- paste0("Cov", seq_len(ncol(covs)))
  }
  if (is.null(rownames(covs))) {
    warning("Missing (row)names in 'covs'.\n  Consistency with omics row names cannot be evaluated")
  }
  else {
    for (l in seq_len(M)) {
      if (!all(dim.names[[1]][[l]]==rownames(covs))) {
        warning("Row names across omic blocks are inconsistent.")
      }
    }
  }
  if (inherits(covs, "data.frame")) {
    form <- stats::as.formula(paste0(" ~ -1 + ",
                                     paste0(colnames(covs), 
                                            collapse = " + ")))
    Q <- stats::model.matrix(form,data=stats::model.frame(form, 
                                                          data = covs,
                                                          na.action="na.pass"))
  }
  
  
  else  {
    Q <- stats::model.matrix(~ -1 + x,stats::model.frame(~ -1 + x,
                                                         data.frame(x=covs[,1]),
                                                         na.action="na.pass"))
  }
  Q <- scale(Q)
  Q[is.na(Q)] <- 0
  SVD.Q <- svd(Q, nv = 0)
  U.Q <- SVD.Q$u[, SVD.Q$d^2 > 1e-05]
  R <- lapply(seq_len(M), function(r) crossprod.na(U.Q, data.blocks[[r]]))
  L <- lapply(R, function(x) crossprod.na(t(U.Q), x))
  block.class <- rep("matrix", M)
  block.class[vapply(data.blocks, inherits, TRUE, what = "FBM")] <- "FBM"
  data.blocks.adj <- lapply(seq_len(M), function(r) {
    if (inherits(data.blocks[[r]], "FBM")) {
      CC <- bigstatsr::FBM(nrow(data.blocks[[r]]), ncol(data.blocks[[r]]),
                           create_bk = T)
      bigstatsr::big_apply(CC, function(x, ind) {
        x[, ind] <- data.blocks[[r]][, ind] - L[[r]][,
                                                     ind]
        NULL
      }, a.combine = "c", ind = seq_along(dim.names[[r]][[2]]))
    }
    else CC <- data.blocks[[r]] - L[[r]]
    return(CC)
  })
  names(data.blocks.adj) <- names(data.blocks)
  return(data.blocks.adj)
}

