# Function to infer clusters separation using the sillouthe index
# (from Taskesen et al).
dbscan_SH <- function(data, eps = NULL,
                      prop_outliers = 0.1,
                      eps_res = 500, eps_range = NULL) {
  data <- data.frame(data)
  if (is.null(eps)) {
    cat("Optimising eps: \n")
    if (is.null(eps_range)) {
      eps_scale <- mean(apply(data, 2, stats::sd))
      epsvec <- seq(0, 4, length.out = eps_res) * eps_scale
    } else {
      epsvec <- seq(eps_range[1], eps_range[2], length.out = eps_res)
    }
    silvec <- numeric(length(epsvec))
    for (i in seq_len(length(epsvec))) {
      eps <- epsvec[i]
      DBcl <- dbscan::dbscan(data, eps)
      cl <- DBcl$cluster
      cat("  Tuning cluster partition: iteration", i, "of", length(epsvec), "\n")
      if (all(cl == 1)) {
        break
      } else if (max(cl) == 1) {
        silvec[i] <- 0
      } else
      if (all(cl == 0)) {
        silvec[i] <- 0
      } else
      if (mean(cl == 0) > prop_outliers) {
        silvec[i] <- 0
      } else {
        S <- cluster::silhouette(
          x = cl[cl != 0],
          dist = stats::dist(data[cl != 0, ])
        )
        silvec[i] <- summary(S)$avg.width
      }
    }
    cat("\n")
    eps <- epsvec[which.max(silvec)]
  }
  DBcl <- dbscan::dbscan(data, eps)$cluster
  cat("Used eps: ", eps, "\n")
  return(list(cluster = DBcl, eps = eps, SIL = max(silvec)))
}
