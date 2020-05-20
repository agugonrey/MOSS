#' Sparse Singular Value Decomposition via Elastic Net.
#'
#' This function performs sparse SVD at subjects' level, features, or both, via Elastic Net types of penalties.
#' It uses a penalized SVD, or PCA, that imposes sparsity in both, left and right, singular vectors.
#' For one PC (rank 1 case), the algorithm finds vectors u, w that minimize:
#'    ||x - u w'||_F^2 + lambda_w (alpha_w||w||_1 + (1 - alpha_w)||w||_F^2) + lambda_u (alpha||u||_1 + (1 - alpha_u)||u||_F^2)
#' such that ||u|| = 1. The right Eigen vector is obtained from v = w / ||w|| and the corresponding Eigen value = u^T x v.
#' The penalties lambda_u and lambda_w are mapped from specified desired degree of sparsity (dg.spar.features & dg.spar.subjects).
#'
#' The function allows the use of the base svd function for relatively small problems. For larger problems, functions for fast-partial SVD (irlba and big_randomSVD, from irlba and bigstatsr packages)
#' are used.
#'
#' @param O Numeric matrix of n subjects (rows) and p features (columns). It can be a Filebacked Big Matrix.
#' @param n.PC Number of desired principal axes. Numeric. Defaults to 1.
#' @param dg.spar.features Degree of sparsity at the features level. Numeric. Defaults to NULL.
#' @param dg.spar.subjects Degree of sparsity at the subjects level. Numeric. Deafults to NULL.
#' @param maxit Maximum number of iterations for the sparse SVD algorithm. Numeric. Defaults to 500.
#' @param tol Convergence tolerance for the sparse SVD algorithm. Numeric. Defaults to 0.001.
#' @param alpha.f Elastic net mixture parameter at the features level. Measures the compromise between lasso (alpha = 1) and ridge (alpha = 0) types of sparsity. Numeric. Deafaults to 1.
#' @param alpha.s Elastic net mixture parameter at the subjects level. Defaults to alpha.s = 1.
#' @param center.arg Should O be centered? Logical. Defaults to TRUE.
#' @param scale.arg Should O be scaled? Logical. Defaults to TRUE.
#' @param approx.arg Should we use standard SVD or random approximations? Defaults to FALSE. If TRUE & is(O,'matrix') == TRUE, irlba is called. If TRUE & is(O, "FBM") == TRUE, big_randomSVD is called.
#' @param svd.0 List containing an initial SVD. Defaults to NULL.
#' @param s.values Should the singular values be calculated? Logical. Defaults to TRUE.
#' @param ncores Number of cores used by big_randomSVD. Default does not use parallelism. Ignored when class(O)!=FBM.
#' @return A list with the results of the (sparse) SVD, containing:
#' \itemize{
#' \item u: Matrix with left eigenvectors.
#' \item v: Matrix with right eigenvectors.
#' \item d: Matrix with singular values. 
#' }
#' @export
#' @references
#'  \itemize{
#'    \item Shen, Haipeng, and Jianhua Z. Huang. 2008. Sparse Principal Component Analysis via Regularized Low Rank Matrix Approximation. Journal of Multivariate Analysis 99 (6). Academic Press:1015_34.
#'    \item Baglama, Jim, Lothar Reichel, and B W Lewis. 2018. Irlba: Fast Truncated Singular Value Decomposition and Principal Components Analysis for Large Dense and Sparse Matrices.
#'  }
#' @examples
#' library("MOSS")
#'
#' #Extracting simulated omic blocks.
#' sim_blocks <- simulate_data()$sim_blocks
#' 
#' X <- sim_blocks$`Block 3`
#'
#' #Equal to svd solution: exact singular vectors and values.
#' out <- ssvdEN(X,approx.arg = FALSE)
#' 
#' \dontrun{
#' #Uses irlba to get approximated singular vectors and values.
#' library(irlba)
#' out <- ssvdEN(X, approx.arg = TRUE)
#' #Uses bigstatsr to get approximated singular vectors and values of a Filebacked Big Matrix.
#' library(bigstatsr)
#' out <- ssvdEN(as_FBM(X), approx.arg = TRUE)
#' }
#' 
#' #Sampling a number of subjects and features for a fix sparsity degree.
#' s.u <- sample(1:nrow(X), 1)
#' s.v <- sample(1:ncol(X), 1)
#'
#' #Lasso penalties.
#' all.equal(sum(ssvdEN(X,dg.spar.features = s.v)$v != 0),s.v)
#' all.equal(unique(colSums(ssvdEN(X,dg.spar.features = s.v,n.PC=5)$v != 0)), s.v)
#'
#' all.equal(sum(ssvdEN(X,dg.spar.subjects  = s.u)$u != 0),s.u)
#' all.equal(unique(colSums(ssvdEN(X,dg.spar.subjects = s.u,n.PC=5)$u != 0)), s.u)
#'
#' out <- ssvdEN(X,dg.spar.features = s.v,dg.spar.subjects = s.u)
#' all.equal(sum(out$u != 0), s.u)
#' all.equal(sum(out$v != 0), s.v)
#'
#' out <- ssvdEN(X,dg.spar.features = s.v,dg.spar.subjects = s.u,n.PC=10)
#' all.equal(unique(colSums(out$u != 0)), s.u)
#' all.equal(unique(colSums(out$v != 0)), s.v)
#'
#' #Ridge penalties.
#' all.equal(sum(ssvdEN(X,dg.spar.features = s.v,alpha.f = 0)$v != 0),ncol(X))
#' all.equal(unique(colSums(ssvdEN(X,dg.spar.features = s.v, n.PC=5,alpha.f = 0)$v != 0)), ncol(X))
#'
#' all.equal(sum(ssvdEN(X,dg.spar.subjects = s.u,alpha.s = 0)$u != 0),nrow(X))
#' all.equal(unique(colSums(ssvdEN(X,dg.spar.subjects = s.u, n.PC=5,alpha.s = 0)$u != 0)), nrow(X))
#'
#' out <- ssvdEN(X,dg.spar.features = s.v,dg.spar.subjects = s.u,alpha.f = 0,alpha.s = 0)
#' all.equal(sum(out$u != 0), nrow(X))
#' all.equal(sum(out$v != 0), ncol(X))
#'
#' out <- ssvdEN(X,dg.spar.features = s.v,dg.spar.subjects = s.u,n.PC=10,alpha.f = 0,alpha.s = 0)
#' all.equal(unique(colSums(out$u != 0)), nrow(X))
#' all.equal(unique(colSums(out$v != 0)), ncol(X))
#' 
#' #Elastic Net penalties.
#' sum(ssvdEN(X,dg.spar.features = s.v,alpha.f = 0.5)$v != 0.5) >= s.v
#' all(unique(colSums(ssvdEN(X,dg.spar.features = s.v, n.PC=5,alpha.f = 0.5)$v != 0)) >= s.v)
#'
#' sum(ssvdEN(X,dg.spar.subjects = s.u,alpha.s = 0.5)$v != 0) >= s.u
#' all(unique(colSums(ssvdEN(X,dg.spar.subjects = s.u, n.PC=5,alpha.s = 0.5)$u != 0)) >= s.u)
#' 
#' #Example of usage within moss.
#'
#' out <- moss(sim_blocks[-4],
#'      K.X=1,
#'      dg.grid.right = 22,
#'      dg.grid.left = 11,
#'      alpha.right = 1,
#'      alpha.left = 1)
#'
#' colSums(out$sparse$u!=0)
#' colSums(out$sparse$v!=0)
ssvdEN <- function (O, n.PC = 1, dg.spar.features = NULL, dg.spar.subjects = NULL,maxit = 500, tol = 0.001,
                    scale.arg = TRUE, center.arg = TRUE, approx.arg = FALSE, alpha.f = 1, alpha.s = 1, svd.0 = NULL,s.values=TRUE, ncores=1) {

  #Checking if the right packages are present to handle approximated SVDs.
  if (approx.arg == TRUE) {
    if (inherits(O, "FBM") == TRUE)  {if (!requireNamespace("bigstatsr",quietly = TRUE)) stop("Package bigstatsr needs to be installed to handle FBM objects.")}
    else {if (!requireNamespace("irlba",quietly = TRUE)) stop("Package irlba needs to be installed to get fast truncated SVD solutions.")}
  }
  
  #Getting matrix dimensions.
  n <- nrow(O)
  p <- ncol(O)

  #If dg.spar.subjects = NULL, get full subjects' scores.
  if (is.null(dg.spar.subjects) == TRUE) dg.spar.subjects <- n

  #If dg.spar.features = NULL, get full features' scores.
  if (is.null(dg.spar.features) == TRUE) dg.spar.features <- p

  if (is.null(svd.0) == TRUE) {
    if (approx.arg == TRUE) {
      if (inherits(O, "FBM") == TRUE) s <- bigstatsr::big_randomSVD(O, fun.scaling = bigstatsr::big_scale(center = center.arg, scale = scale.arg), k = n.PC, ncores = ncores)
      else {
        O <- scale(O, center=center.arg, scale = scale.arg)
        s <- irlba::irlba(O, nu = n.PC, nv = n.PC)[c("u","v","d")]
      }
    }
    else {
      O <- scale(O, center=center.arg, scale = scale.arg)
      s <- svd(O, nu = n.PC, nv = n.PC); s$d <- s$d[1 : n.PC]
    }
  }
  else s <- svd.0

  if (dg.spar.features <= p || dg.spar.subjects <= n) {
    dg.spar.features <- p - dg.spar.features
    dg.spar.subjects <- n - dg.spar.subjects

    if (length(dg.spar.features) != n.PC) dg.spar.features <- rep(dg.spar.features, length.out = n.PC)
    if (length(dg.spar.subjects) != n.PC) dg.spar.subjects <- rep(dg.spar.subjects, length.out = n.PC)

    #Elastic net solutions.
    softEN <- function(y, dg.spar, alpha) {
      a <- abs(y)
      z <- sort(a)

      #Mapping dg to lambda scale.arg
      lambda <- z[dg.spar]
      b <- a - lambda * alpha
      y.out <- sign(y) * ifelse(b > 0, b, 0) / (1 + (1 - alpha) * lambda)
      
      return(y.out)
    }
    for(j in 1:n.PC) s$v[,j] <- s$d[j] * s$v[, j]

    #Coordinate descend: iterates until convergence or until reaching maxit.
    iter <- 1
    delta_u <- Inf
    shrink.features <- 1 : length(dg.spar.features)
    shrink.features[which(dg.spar.features == 0)] <- 0
    shrink.subjects <- 1 : length(dg.spar.subjects)
    shrink.subjects[which(dg.spar.subjects == 0)] <- 0

    while (delta_u > tol && iter < maxit) {
      u <- s$u
      usx <- t(s$u) %*% O
      s$v <- qr.Q(qr(t(usx)))

      #Solving v for fixed u.
      if (sum(shrink.features) > 0) for (j in shrink.features) s$v[, j] <- softEN(s$v[, j], dg.spar.features[j], alpha.f)
      norm. <- stats::sd(s$v)
      if (norm. == 0) norm. <- 1
      s$v <- s$v / norm.

      #Getting orthogonal subject's scores.
      xsv <- O %*% s$v
      s$u <- qr.Q(qr(xsv))

      #Solving u for fixed v.
      if (sum(shrink.subjects) > 0) for (j in shrink.subjects) s$u[, j] <- softEN(s$u[, j], dg.spar.subjects[j], alpha.s)
      norm. <- stats::sd(s$u)
      if (norm. == 0) norm. <- 1
      s$u <- s$u / norm.

      #Fixing direction
      s$u <- sweep(s$u, 2, 
                   apply(xsv, 2, function(x) sign(utils::head(O[O[,1] != 0,1], 1)))/
                     apply(s$u, 2, function(x) sign(utils::head(O[O[,1] != 0,1], 1))), 
                   `*`)
      
      #Calculating convergence rate in terms of Frobenious norm.
      delta_u <- norm(u - s$u, "f") / (n * n.PC)
      iter <- iter + 1
    }
    if (iter >= maxit) warning("Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.")

    #Getting scales ("Eigen" values).
    s$v <- tcrossprod(s$v, diag(1/sqrt(apply(s$v, 2, crossprod)), ncol(s$v), ncol(s$v)))
    if (s.values == TRUE) s$d <- crossprod(s$u, O %*% s$v)
    else s$d <- NULL
  }
  return(s)
}
