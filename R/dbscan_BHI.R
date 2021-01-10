# Function to infer clusters separation by
# maximizing the Biological Homogeneity Index.
dbscan_BHI <- function(data,
                       eps = NULL,
                       minPts = 5,
                       prop_outliers = 0.1,
                       eps_res = 500,
                       eps_range = NULL,
                       annot, plot = F) {
  data <- data.frame(data)
  if (is.null(eps)) {
    cat("Optimising eps: ")
    if (is.null(eps_range)) {
      eps_scale <- mean(apply(data, 2, stats::sd))
      # makes the search scale independent
      epsvec <- seq(0, 4, length.out = eps_res) * eps_scale
      # space to search for eps parameter
    } else {
      epsvec <- seq(eps_range[1], eps_range[2], length.out = eps_res)
    }
    bhi <- numeric(length(epsvec))
    if (plot) {
      graphics::plot(
        x = epsvec, y = seq(0, 1,
          length.out = length(epsvec)
        ),
        pch = NA, xlab = "EPS", ylab = "BHI"
      )
    }
    for (i in seq_len(length(epsvec))) {
      eps <- epsvec[i]
      DBcl <- dbscan::dbscan(data, eps, MinPts = minPts) # quite fast
      cl <- DBcl$cluster
      cat(".")
      if (all(cl == 1)) {
        break
      } else if (max(cl) == 1) {
        bhi[i] <- 0
      } else
      if (all(cl == 0)) {
        bhi[i] <- 0
      } else
      if (mean(cl == 0) > prop_outliers) {
        bhi[i] <- 0
      } else {
        bhi[i] <- BHI(
          statClust = cl[cl != 0],
          annotation = annot[cl != 0, ]
        ) # exclude the 0'
      }
      if (plot) graphics::points(y = bhi[i], x = eps, pch = 16, col = 2)
    }
    eps <- epsvec[which.max(bhi)]
    bhi <- bhi[which.max(bhi)]
    cat("\n")
  }
  DBcl <- dbscan::dbscan(data, eps, MinPts = minPts)$cluster

  if (!is.null(eps)) {
    bhi <- BHI(
      statClust = DBcl[DBcl != 0],
      annotation = annot[DBcl != 0, ]
    )
  } # exclude the 0'

  cat("Used eps: ", eps, "( BHI = ", round(bhi, 2), ")\n")
  return(list(cluster = DBcl, eps = eps, BHI = bhi))
}
