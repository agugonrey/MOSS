#' 'Solution path' for sparse Singular Value Decomposition via Elastic Net
#' using parallel computing.
#'
#' This function is a copy of 'ssvdEN_sol_path' meant to be used
#' in combination with the future.apply package to allow for
#' parallel computing of the optimal degrees of sparsity by
#' subjects and/or features.
#' 
#' @note Although the degree of sparsity maps onto number of
#' features/subjects for Lasso, the user needs to be aware that
#' this conceptual correspondence
#' is lost for full EN (alpha belonging to (0, 1);
#' e.g. the number of features selected with alpha < 1
#' will be eventually larger than the optimal degree of sparsity).
#' This allows to rapidly increase the number of non-zero elements
#' when tuning the degrees of sparsity.
#' In order to get exact values for the degrees of sparsity at subjects or
#'  features levels, the user needs to
#' set the value of 'exact.dg' parameter from 'FALSE' (the default) to
#' 'TRUE'.
#' @param O Numeric matrix of n subjects (rows) and p features (columns).
#'  Only objects supported are 'matrix' and 'FBM'.
#' @param n.PC Number of desired principal axes. Numeric. Defaults to 1.
#' @param dg.grid.right Grid with degrees of sparsity at the features level.
#' Numeric. Default is the entire solution path for features
#' (i.e. 1 : (ncol(O) - 1)).
#' @param dg.grid.left Grid with degrees of sparsity at the subjects level.
#' Numeric. Defaults to dg.grid.left = nrow(O).
#' @param svd.0 Initial SVD (i.e. least squares solution).
#' Defaults to NULL.
#' @param alpha.f Elastic net mixture parameter at the features level.
#' Measures the compromise between lasso (alpha = 1) and
#'  ridge (alpha = 0) types of sparsity. Numeric. Defaults to 1.
#' @param alpha.s Elastic net mixture parameter at the subjects level.
#'  Defaults to alpha.s = 1.
#' @param tol Convergence is determined when ||U_j - U_{j-1}||_F < tol,
#'  where U_j is the matrix of estimated left regularized singular
#'  vectors at iteration j.
#' @param center Should we center? Logical. Defaults to TRUE.
#' @param scale Should we scale? Logical. Defaults to TRUE.
#' @param maxit Maximum number of iterations. Defaults to 500.
#' @param ncores Number of cores used by big_randomSVD.
#' Default does not use parallelism. Ignored when is(O, "FBM") == TRUE.
#' @param plot Should we plot the solution path? Logical. Defaults to FALSE
#' @param right.lab Label for the features level. Character.
#' Defaults to 'features'.
#' @param left.lab Label for the subjects level. Character.
#' Defaults to 'subjects'.
#' @param approx Should we use standard SVD or random approximations?
#' Defaults to FALSE. If TRUE & is(O,'matrix') == TRUE, irlba is called.
#'  If TRUE & is(O, "FBM") == TRUE, big_randomSVD is called.
#' @param verbose Should we print messages?. Logical. Defaults to TRUE.
#' @param exact.dg Should we compute exact degrees of sparsity? Logical.
#' Defaults to FALSE. Only relevant When alpha.s or alpha.f are in the (0,1)
#'  interval and exact.dg = TRUE.
#' @param lib.thresh Should we use a liberal or conservative
#' threshold to tune degrees of sparsity? Logical. Defaults to TRUE.
#' @export
#' @examples
#'\donttest{
#' library("MOSS")
#'
#' # Extracting simulated omic blocks.
#' sim_blocks <- simulate_data()$sim_blocks
#' X <- sim_blocks$`Block 3`
#'
#' # Comparing ssvdEN_sol_path_par and ssvdEN_sol_path.
#' t1 <- proc.time()
#' out1 <- ssvdEN_sol_path(X, dg.grid.right = 1:1000, dg.grid.left = 1:500)
#' t1 <- proc.time() - t1
#' 
#' t2 <- proc.time()
#' out2 <- ssvdEN_sol_path_par(X, dg.grid.right = 1:1000, dg.grid.left = 1:500)
#' t2 <- proc.time() - t2
#' }
ssvdEN_sol_path_par <- function(O, center = TRUE, scale = TRUE,
                            dg.grid.right = seq_len(ncol(O)) - 1,
                            dg.grid.left = NULL,
                            n.PC = 1, svd.0 = NULL,
                            alpha.f = 1, alpha.s = 1,
                            maxit = 500, tol = 1E-03,
                            approx = FALSE, plot = FALSE, ncores = 1,
                            verbose = TRUE,lib.thresh=TRUE,
                            left.lab = "Subjects",
                            right.lab = "Features",
                            exact.dg = FALSE) {

  # Checking if the right packages.
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package future.apply needs to be installed to perform parallel
         computing.")
  }
  
  if (inherits(O, "FBM") == TRUE) {
    if (!requireNamespace("bigstatsr", quietly = TRUE)) {
      stop("Package bigstatsr needs to be installed to handle FBM objects.")
    }
    if (approx == FALSE) approx <- TRUE
  }
  else {
      if (!requireNamespace("irlba", quietly = TRUE)) {
        stop("Package irlba needs to be installed to get 
           fast truncated SVD solutions.")
      }
  }

  if (any(vapply(c("matrix", "array", "FBM"),
    FUN =
      function(x) inherits(O, x), TRUE
  )) == FALSE) {
    stop("Input needs to be a matrix or FBM.")
  }
  if (!(any(dg.grid.right %in% seq_len(ncol(O))) |
    is.null(dg.grid.right)) |
    !(any(dg.grid.left %in% seq_len(nrow(O))) | is.null(dg.grid.left))) {
    stop(
      "Degrees of sparsity for ", right.lab, " and ", left.lab,
      " need to be NULL or belong to [1, ", ncol(O), "]",
      " and ", "[1, ", nrow(O), "], respectively."
    )
  }

  # Checking if the right packages are present for plotting.
  if (plot == TRUE) {
    if (!requireNamespace("viridis", quietly = TRUE)) {
      stop("Package 'viridis' needs to be installed for 
             graphical displays.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' needs to be installed for
             graphical displays.")
    }
    if (!requireNamespace("ggpmisc", quietly = TRUE)) {
      stop("Package 'ggpmisc' needs to be installed for
             showing peaks on the PEV trajectory.")
    }
  }

  if (is.null(svd.0) == TRUE) {
    if (approx == TRUE) {
      if (inherits(O, "FBM") == TRUE) {
        s <- bigstatsr::big_randomSVD(O,
          fun.scaling = bigstatsr::big_scale(
            center = center,
            scale = scale
          ),
          k = n.PC, ncores = ncores
        )
      } 
      else {
        O <- scale(O, center = center, scale = scale)
        s <- irlba::irlba(O, nu = n.PC, nv = n.PC)[c("u", "v", "d")]
      }
    }
    else {
      O <- scale(O, center = center, scale = scale)
      s <- svd(O, nu = n.PC, nv = n.PC)
      s$d <- s$d[1:n.PC]
    }
  }
  else {
    s <- svd.0
    s$u <- as.matrix(s$u)
    s$v <- as.matrix(s$v)
  }

  # At what levels shall we tune?
  n.s <- length(dg.grid.left)
  n.f <- length(dg.grid.right)

  # Tuning parameters.
  PEV <- NULL
  if (n.s > 1) {
    # Getting proportion of explained variance by degree of sparsity.
    pev <- do.call("c", future.apply::future_lapply(1:n.s, function(i) {
      if (verbose) {
        cat(paste0(
          "Tuning ",
          left.lab,
          " degree of sparsity = ",
          dg.grid.left[i], "  (max value on the grid ",
          dg.grid.left[n.s], ").\n"
        ))
      }
      s1 <- ssvdEN(O,
        svd.0 = s,
        n.PC = n.PC,
        dg.spar.subjects = dg.grid.left[i],
        dg.spar.features = NULL,
        alpha.s = alpha.s,
        tol = tol,
        maxit = maxit,
        ncores = ncores,
        exact.dg = exact.dg
      )
      sum(diag(crossprod(s1$u %*% s1$d)))
    }))
    pev <- pev / max(pev)
    # Getting first and second empirical derivatives.
    dpev1 <- c(0, diff(pev) / diff(dg.grid.left))
    dpev2 <- c(0, diff(dpev1) / diff(dg.grid.left))
    PEV$u <- data.frame(
      y = c(pev, dpev1, -dpev2),
      x = rep(dg.grid.left, times = 3),
      type = rep(c(
        "PEV",
        "PEV first derivative",
        "PEV second reciprocal derivative"
      ),
      each = length(pev)
      ),
      facet = rep(left.lab, length(pev))
    )

    if (lib.thresh) opt.dg.left <- dg.grid.left[which.max(dpev1)]
    else {
      th <- which(dpev2 ^ 2 >= mean(dpev2 ^ 2) + stats::sd(dpev2 ^ 2))
      nf_dsvd_d2 <- ifelse(th[length(th)] + 2 <= max(dg.grid.left),
                           th[length(th)] + 2, 
                           max(dg.grid.left))
      opt.dg.left <- dg.grid.left[nf_dsvd_d2]
    }
  }
  # If a single or no value for dg of spar for subjects is given...
  else {
    opt.dg.left <- dg.grid.left
  }

  if (n.f > 1) {
    # Getting proportion of explained variance by degree of sparsity.
    pev <- do.call("c", future.apply::future_lapply(1:n.f, function(j) {
      if (verbose) {
        cat(paste0(
          "Tuning ", right.lab,
          " degree of sparsity = ",
          dg.grid.right[j],
          "  (max value on the grid= ",
          dg.grid.right[n.f], ").\n"
        ))
      }

      s1 <- ssvdEN(O,
        svd.0 = s,
        n.PC = n.PC,
        dg.spar.subjects = NULL,
        dg.spar.features = dg.grid.right[j],
        alpha.f = alpha.f,
        tol = tol,
        maxit = maxit,
        ncores = ncores,
        exact.dg = exact.dg
      )
      sum(diag(crossprod(s1$u %*% s1$d)))
    }))
    pev <- pev / max(pev)

    # Getting first and second empirical derivatives.
    dpev1 <- c(0, diff(pev) / diff(dg.grid.right))
    dpev2 <- c(0, diff(dpev1) / diff(dg.grid.right))
    PEV$v <- data.frame(
      y = c(pev, dpev1, -dpev2),
      x = rep(dg.grid.right, times = 3),
      type = rep(c(
        "PEV",
        "PEV first derivative",
        "PEV second reciprocal derivative"
      ),
      each = length(pev)
      ),
      facet = rep(right.lab, length(pev))
    )

    if (lib.thresh) opt.dg.right <- dg.grid.right[which.max(dpev1)]
    else {
      th <- which(dpev2 ^ 2 >= mean(dpev2 ^ 2) + stats::sd(dpev2 ^ 2))
      nf_dsvd_d2 <- ifelse(th[length(th)] + 2 <= max(dg.grid.right),
                           th[length(th)] + 2, 
                           max(dg.grid.right))
      opt.dg.right <- dg.grid.right[nf_dsvd_d2]
    }
    
  }
  # If a single or no value for dg of spar for features is given...
  else {
    opt.dg.right <- dg.grid.right
  }

  # Preparing outcome.
  out <- NULL

  if ((n.f > 1 || n.s > 1) && plot) {
    # Creating a ggplot object.
    # Preparing plots.
    d <- do.call("rbind", PEV)
    suppressWarnings(out$plot <- ggplot2::ggplot(
      d,
      ggplot2::aes_string(x = "x", y = "y")
    ) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(facet ~ type, scales = "free") +
      ggpmisc::stat_peaks(colour = "#FCA50AB3") +
      ggpmisc::stat_peaks(
        geom = "text",
        colour = "#2A788EB3", angle = 30,
        hjust = -0.5
      ) +
      ggplot2::scale_x_continuous("Degrees of sparsity") +
      ggplot2::scale_y_continuous("Proportion of explained variance 
                                   \nand derivatives") +
      ggplot2::theme_minimal())
  }
  else {
    out$plot <- NULL
  }


  # Returning SVD for the selected features.
  out$SVD <- ssvdEN(O,
    svd.0 = s,
    n.PC = n.PC,
    dg.spar.subjects = opt.dg.left,
    dg.spar.features = opt.dg.right,
    alpha.s = alpha.s,
    alpha.f = alpha.f,
    tol = tol,
    exact.dg = exact.dg,
    maxit = maxit,
    ncores = ncores
  )

  # Adding the solutions into the output.
  out$SVD$opt.dg.left <- opt.dg.left
  out$SVD$opt.dg.right <- opt.dg.right
  return(out)
}
