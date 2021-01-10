# PCA of an extended data matrix.
rrr_pca <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {
  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))

  # If any data block is of class "FBM", the extended data matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p)
    b.names <- rep(names(data.blocks), p)
    for (m in 1:M) {
      if (verbose) message(" Block ", m, " of ", M)
      # Within blocks chunks.
      block.chunks <-
        rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
          each = bigstatsr::block_size(p[m], 1)
        )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names ==
          names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks)
  }
  s <- rrr_big(X = Z, K.X = K.X, verbose = verbose, ncores = ncores)
  return(list("SVD" = s, "Z" = Z))
}

# Multi-block PCA of a series of data blocks.
rrr_mb_pca <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {

  # If M = 1, do conventional PCA.
  if (M == 1) {
    return(rrr_pca(data.blocks))
  }

  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks), p)

  # If any data block is of class "FBM",
  # the extended data matrix is turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p)
    for (m in 1:M) {
      if (verbose) message(" Block ", m, " of ", M)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }

  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks)
  }

  # Creating a list to store the SVD results.
  SVD <- vector("list", M + 1)
  SVD[[1]]$d <- 0
  SVD[[1]]$v <- matrix(0, sum(p), K.X)

  # Calculating SVD by block.
  for (m in 2:(M + 1)) {
    if (verbose) message("SVD for block ", m - 1, " of ", M)
    SVD[[m]] <- rrr_big(
      X = data.blocks[[m - 1]], K.X = K.X,
      verbose = verbose, ncores = ncores
    )
    SVD[[1]]$d <- SVD[[1]]$d + SVD[[m]]$d
    SVD[[1]]$v[which(b.names == names(data.blocks)[m - 1]), ] <- SVD[[m]]$v
  }

  # Calculating block weights.
  for (m in 2:(M + 1)) SVD[[m]]$b <- SVD[[m]]$d / SVD[[1]]$d

  # Getting global matrix L.Z.
  if (verbose) message("SVD for extended data matrix.")
  if (inherits(Z, "FBM") == TRUE) {
    count <- 0
    SVD[[1]]$u <-
      bigstatsr::big_apply(Z, function(X, ind) {
        if (verbose) {
          cat(
            " Working chunk ", (count <<- count + 1), "of",
            ceiling(ncol(Z) / bigstatsr::block_size(ncol(Z), 1)),
            ".\n"
          )
        }
        Z[, ind] %*% SVD[[1]]$v[ind, ]
      },
      ind = seq_len(ncol(Z)),
      a.combine = "plus",
      block.size = bigstatsr::block_size(ncol(Z), 1)
      )
  }
  else {
    SVD[[1]]$u <- Z %*% SVD[[1]]$v
  }
  # Rotating u.
  if (verbose) message("Generating orthogonal axes.")
  SVD[[1]]$u <- qr.Q(qr(SVD[[1]]$u))
  return(list("SVD" = SVD, "Z" = Z))
}

# PLS of a multivariate response and an extended matrix of predictors.
rrr_pls <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {
  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }
  s <- rrr_big(
    X = Z, Y = data.blocks[[1]], K.X = K.X, K.Y = K.Y,
    verbose = verbose, ncores = ncores
  )
  return(s)
}

# Multi-block PLS of a series of data blocks.
rrr_mb_pls <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {

  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix
  # is turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }

  # Creating a list to store the SVD results.
  SVD <- vector("list", M)
  SVD[[1]]$d <- 0
  SVD[[1]]$v <- matrix(0, ncol(Z), K.Y)

  # Calculating SVD by block.
  for (m in 2:M) {
    if (verbose) message("PLS ", m - 1, " of ", M - 1)
    SVD[[m]] <- rrr_big(
      Y = data.blocks[[1]],
      X = data.blocks[[m]],
      K.X = K.X, K.Y = K.Y,
      verbose = verbose, ncores = ncores
    )$SVD
    SVD[[1]]$d <- SVD[[1]]$d + SVD[[m]]$d
    SVD[[1]]$v[which(b.names == names(data.blocks)[m]), ] <- SVD[[m]]$u
  }

  # Calculating block weights.
  for (m in 2:M) SVD[[m]]$b <- SVD[[m]]$d / SVD[[1]]$d

  if (verbose) message("Getting global LR matrix.")
  LR <- rrr_big(
    Y = Z, X = data.blocks[[1]],
    K.X = K.Y, K.Y = K.X,
    verbose = verbose,
    ncores = ncores,
    lr.return = T
  )
  if (verbose) message("Getting latent factors for response block.")

  if (inherits(LR, "FBM") == TRUE) {
    count <- 0
    SVD[[1]]$u <-
      bigstatsr::big_apply(LR, function(X, ind) {
        if (verbose) {
          cat(
            "Working chunk ",
            (count <<- count + 1), "of",
            ceiling(ncol(LR) / bigstatsr::block_size(
              ncol(LR),
              1
            )), ".\n"
          )
        }
        LR[, ind] %*% SVD[[1]]$v[ind, ]
      },
      ind = seq_len(ncol(LR)),
      a.combine = "plus",
      block.size = bigstatsr::block_size(ncol(LR), 1)
      )
  }
  else {
    SVD[[1]]$u <- LR %*% SVD[[1]]$v
  }
  # Rotating u.
  if (verbose) message("Generating orthogonal axes.")
  SVD[[1]]$u <- qr.Q(qr(SVD[[1]]$u))
  aux <- SVD[[1]]
  SVD[[1]]$v <- aux$u
  SVD[[1]]$u <- aux$v
  return(list("SVD" = SVD, "LR" = LR))
}

# RRR of a multivariate response and an extended matrix of predictors.
rrr_rrr <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {
  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }
  s <- rrr_big(
    X = Z, Y = data.blocks[[1]], power = -1, K.X = K.X, K.Y = K.Y,
    verbose = verbose, ncores = ncores
  )
  return(s)
}

# Multi-block PLS of a series of data blocks.
rrr_mb_rrr <- function(data.blocks, block.class, K.X, K.Y, ncores, verbose, M) {

  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }

  # Creating a list to store the SVD results.
  SVD <- vector("list", M)
  SVD[[1]]$d <- 0
  SVD[[1]]$v <- matrix(0, ncol(Z), K.Y)

  # Calculating SVD by block.
  for (m in 2:M) {
    if (verbose) message("PLS ", m - 1, " of ", M - 1)
    SVD[[m]] <- rrr_big(
      Y = data.blocks[[1]],
      X = data.blocks[[m]],
      K.X = K.X, K.Y = K.Y,
      verbose = verbose, ncores = ncores,
      power = -1
    )$SVD
    SVD[[1]]$d <- SVD[[1]]$d + SVD[[m]]$d
    SVD[[1]]$v[which(b.names == names(data.blocks)[m]), ] <- SVD[[m]]$u
  }

  # Calculating block weights.
  for (m in 2:M) SVD[[m]]$b <- SVD[[m]]$d / SVD[[1]]$d

  if (verbose) message("Getting global LR matrix.")
  LR <- rrr_big(
    Y = Z, X = data.blocks[[1]],
    K.X = K.Y, K.Y = K.X,
    verbose = verbose,
    ncores = ncores,
    power = -1,
    lr.return = T
  )
  if (verbose) message("Getting latent factors for response block.")

  if (inherits(LR, "FBM") == TRUE) {
    count <- 0
    SVD[[1]]$u <-
      bigstatsr::big_apply(LR, function(X, ind) {
        if (verbose) {
          cat(
            "Working chunk ",
            (count <<- count + 1), "of",
            ceiling(ncol(LR) / bigstatsr::block_size(
              ncol(LR),
              1
            )), ".\n"
          )
        }
        LR[, ind] %*% SVD[[1]]$v[ind, ]
      },
      ind = seq_len(ncol(LR)),
      a.combine = "plus",
      block.size = bigstatsr::block_size(ncol(LR), 1)
      )
  }
  else {
    SVD[[1]]$u <- LR %*% SVD[[1]]$v
  }
  # Rotating u.
  if (verbose) message("Generating orthogonal axes.")
  SVD[[1]]$u <- qr.Q(qr(SVD[[1]]$u))
  return(list("SVD" = SVD, "LR" = LR))
}

# LDA of an extended data matrix.
rrr_lda <- function(data.blocks, block.class, K.X, ncores, verbose, M) {

  # Checking that the response block contains binary elements.
  data.blocks[[1]] <- apply(data.blocks[[1]], 2, function(x) {
    x.elements <- table(x)
    n.x.elements <- sort(x.elements)
    aux <- as.numeric(as.character(factor(x,
      levels = names(n.x.elements),
      labels = c(
        -n.x.elements[1],
        n.x.elements[2]
      ) /
        sum(n.x.elements)
    )))
    return(aux)
  })

  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }
  s <- rrr_big(X = Z, K.X = K.X, verbose = verbose, ncores = ncores)
  u <- s$u
  s$w <- crossprod(u, data.blocks[[1]])
  for (k in 1:K.X) s$w[k, ] <- s$w[k, ] / s$d[k]
  s$`w.x` <- s$u %*% s$w
  s$c <- diag(crossprod(crossprod(diag(s$d[1:K.X]), s$w), s$w) / (2 * n))
  s$c2 <- diag(crossprod(s$w, s$w)) / (2 * n)
  return(list("SVD" = s, "Z" = Z))
}

# Multi-block LDA of a series of data blocks.
rrr_mb_lda <- function(data.blocks, block.class, K.X, ncores, verbose, M) {

  # Getting dimensions for each data block.
  n <- nrow(data.blocks[[1]])
  p <- unlist(lapply(data.blocks, ncol))
  b.names <- rep(names(data.blocks)[-1], p[-1])

  # If any data block is of class "FBM", the extended predictors matrix is
  # turned into an "FBM".
  if (any(block.class == "FBM")) {
    # Creating an extended omic matrix.
    if (verbose) message("Creating an extended data matrix.")
    Z <- bigstatsr::FBM(n, sum(p[-1]), create_bk = T)$save()
    ext.mat.slides <- 1:sum(p[-1])
    for (m in 2:M) {
      if (verbose) message(" Predictor block ", m - 1, " of ", M - 1)
      # Within blocks chunks.
      block.chunks <- rep(1:ceiling(p[m] / bigstatsr::block_size(p[m], 1)),
        each = bigstatsr::block_size(p[m], 1)
      )[1:p[m]]
      n.chunks <- length(unique(block.chunks))
      for (chunk in unique(block.chunks)) {
        if (verbose) cat("   Working chunk ", chunk, " of ", n.chunks, "\n")
        Z[, ext.mat.slides[b.names == names(data.blocks)[m]][block.chunks ==
          chunk]] <-
          data.blocks[[m]][, block.chunks == chunk]
      }
    }
  }
  # If all blocks are of class matrix, simply appends data blocks by column.
  else {
    Z <- do.call("cbind", data.blocks[-1])
  }

  # Creating a list to store the SVD results.
  SVD <- vector("list", M)
  SVD[[1]]$d <- 0
  SVD[[1]]$v <- matrix(0, ncol(Z), K.X)

  # Calculating SVD by block.
  for (m in 2:M) {
    if (verbose) message("LDA ", m - 1, " of ", M - 1)
    SVD[[m]] <- rrr_lda(
      data.blocks = data.blocks[c(1, m)],
      block.class = block.class[c(1, m)],
      K.X = K.X,
      verbose = verbose, ncores = ncores, M = 2
    )$SVD
    SVD[[1]]$d <- SVD[[1]]$d + SVD[[m]]$d
    SVD[[1]]$v[which(b.names == names(data.blocks)[m]), ] <- SVD[[m]]$v
  }

  # Calculating block weights and lda function value.
  for (m in 2:M) {
    SVD[[m]]$b <- SVD[[m]]$d / SVD[[1]]$d
    SVD[[m]]$w <- crossprod(diag(SVD[[m]]$b[1:K.X]), SVD[[m]]$w)
  }

  # Getting global matrix L.Z.
  if (verbose) message("SVD for extended data matrix.")
  if (inherits(Z, "FBM") == TRUE) {
    count <- 0
    SVD[[1]]$u <-
      bigstatsr::big_apply(Z, function(X, ind) {
        if (verbose) {
          cat(
            " Working chunk ", (count <<- count + 1), "of",
            ceiling(ncol(Z) / bigstatsr::block_size(ncol(Z), 1)),
            ".\n"
          )
        }
        Z[, ind] %*% SVD[[1]]$v[ind, ]
      },
      ind = seq_len(ncol(Z)),
      a.combine = "plus",
      block.size = bigstatsr::block_size(ncol(Z), 1)
      )
  }
  else {
    SVD[[1]]$u <- Z %*% SVD[[1]]$v
  }
  # Rotating u.
  if (verbose) message("Generating orthogonal axes.")
  SVD[[1]]$u <- qr.Q(qr(SVD[[1]]$u))
  SVD[[1]]$w <- Reduce("+", lapply(SVD[-1], function(x) x$w))
  SVD[[1]]$`w.x` <- SVD[[1]]$u %*% SVD[[1]]$w

  names(SVD) <- c("Global", names(data.blocks)[-1])
  return(list("SVD" = SVD, "Z" = Z))
}
