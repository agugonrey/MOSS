#' Mapping principal components onto a 2D map via tSNE.
#'
#' This function is called by moss whenever 'moss(tSNE=TRUE)' to project
#' latent factors onto two dimensions via t-stochastic neighbor embedding
#'  (tSNE)
#'
#' However, it can be used on any generic data matrix.
#' The function uses the Barnes-Hut tSNE algorithm from Rtsne package,
#'  and uses an iterative procedure to select a tSNE map minimizing the
#'   projection cost across several random initial conditions.
#' The function is inspired by the iterative procedure discussed in
#'  Taskesen et al. 2016 and code originally provided with the publication.
#' @param Z A matrix with axes of variation (typically PCs) as columns
#'  and subjects as rows.
#' @param perp Perplexity value. Defaults to 50.
#' @param n.samples Number of times the algorithm from different random
#'  initial conditions. Defaults to 1.
#' @param n.iter Number of iterations for each run of the algorithm.
#' Defaults to 1000.
#' @return Returns the tSNE (output of function 'Rtsne::Rtsne') with
#' the smallest error across the 'n.samples' random starts.
#' @references \itemize{
#' \item van der Maaten L, Hinton G. Visualizing Data using t-SNE.
#'  J Mach Learn Res. 2008;9: 2579–2605
#'  \item Krijthe JH. Rtsne: T-Distributed Stochastic Neighbor
#'  Embedding using a Barnes-Hut Implementation. 2015
#'  \item Taskesen, Erdogan, Sjoerd M. H. Huisman, Ahmed Mahfouz,
#'     Jesse H. Krijthe, Jeroen de Ridder, Anja van de Stolpe,
#'     Erik van den Akker, Wim Verheagh, and Marcel J. T. Reinders. 2016.
#'      Pan-Cancer Subtyping in a 2D-Map Shows Substructures
#'      That Are Driven by Specific Combinations of Molecular
#'       Characteristics. Scientific Reports 6 (1):24949.
#'
#'  }
#'
#' @export
#' @examples
#' \dontrun{
#' library("MOSS")
#' sim_blocks <- simulate_data()$sim_blocks
#'
#' # Example of pca2tsne usage.
#' Z <- pca2tsne(sim_blocks$`Block 3`, perp = 50, n.samples = 1, n.iter = 1e3)$Y
#' plot(Z, xlab = "x_tSNE(X)", ylab = "y_tSNE(X)")
#'
#' # Example of use within moss.
#' moss(sim_blocks[-4],
#'   tSNE = list(
#'     perp = 50,
#'     n.samples = 1,
#'     n.iter = 1e3
#'   ),
#'   plot = TRUE
#' )$tSNE_plot
#' }
pca2tsne <- function(Z, perp = 50, n.samples = 1, n.iter = 1000) {
  if (is.null(rownames(Z))) {
    sample.names <- paste0(seq_len(nrow(Z)))
  } else {
    sample.names <- rownames(Z)
  }

  message("Embedding ", ncol(Z), " axes onto two dimensions.\n")
  tSNE.out <- vector(mode = "list", length = n.samples)
  cost <- numeric(length = n.samples)
  for (i in seq_len(n.samples)) {
    cat("tSNE sample =", i, "of", n.samples, ".\n")
    tSNE.out[[i]] <- Rtsne::Rtsne(Z,
      dims = 2, perplexity = perp,
      max_iter = n.iter,
      verbose = T,
      pca = F,
      check_duplicates = F
    )
    cost[i] <- tSNE.out[[i]]$costs[which.min(tSNE.out[[i]]$costs)]
  }

  # Getting optimal embedding.
  S.tsne <- tSNE.out[[which.min(cost)]]
  rownames(S.tsne$Y) <- sample.names
  return(S.tsne)
}
