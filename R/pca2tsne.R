#' Mapping principal components onto a 2D map via tSNE.
#'
#' This function is called by moss whenever 'moss(tSNE=TRUE)' to project
#' latent factors onto two dimensions via t-stochastic neighbor embedding
#'  (tSNE) However, it can be used on any generic data matrix.
#' The function uses the Barnes-Hut tSNE algorithm from Rtsne package,
#'  and uses an iterative procedure to select a tSNE map minimizing the
#'   projection cost across several random initial conditions.
#' The function is inspired by the iterative procedure discussed in
#'  Taskesen et al. 2016 and code originally provided with the publication.
#' @param Z A matrix with axes of variation (typically PCs) as columns
#'  and subjects as rows.
#' @param perp Perplexity value. Defaults to 50.
#' @param n.samples Number of times the algorithm starts from 
#' different random initial conditions. Defaults to 1.
#' @param n.iter Number of iterations for each run of the algorithm.
#' @param parallel Should random starts be done in parallel? Logical. 
#' Default to FALSE.
#' Defaults to 1000.
#' @return Returns output of function 'Rtsne::Rtsne' from the random initial
#' condition with the smallest 'reconstruction error'.
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
#' 
#' library("MOSS")
#' sim_blocks <- simulate_data()$sim_blocks
#'
#' # Example of pca2tsne usage.
#' Z <- pca2tsne(sim_blocks$`Block 3`, 
#'               perp = 50, 
#'               n.samples = 1,
#'               n.iter = 1e3)$Y
#' plot(Z, xlab = "x_tSNE(X)", ylab = "y_tSNE(X)")
#'\donttest{
#' # Example of usage within moss.
#' set.seed(34)
#' moss(sim_blocks[-4],
#'   tSNE = list(
#'     perp = 50,
#'     n.samples = 1,
#'     n.iter = 1e3
#'   ),
#'   plot = TRUE
#' )$tSNE_plot
#' }
pca2tsne <- function(Z,
                     perp = 50, 
                     n.samples = 1, 
                     n.iter = 1000, 
                     parallel=FALSE) {
  if (is.null(rownames(Z))) {
    sample.names <- paste0(seq_len(nrow(Z)))
  } else {
    sample.names <- rownames(Z)
  }

  message("Embedding ", ncol(Z), " axes onto two dimensions.\n")
  
  if (n.samples == 1 & parallel == TRUE) parallel <- FALSE
  
  if (parallel == FALSE) {
    tSNE.out <- vector(mode = "list", length = n.samples)
    cost <- numeric(length = n.samples)
    for (i in seq_len(n.samples)) {
      cat("tSNE sample =", i, "of", n.samples, ".\n")
      tSNE.out[[i]] <- Rtsne::Rtsne(Z,
                                    dims = 2, perplexity = perp,
                                    max_iter = n.iter,
                                    verbose = TRUE,
                                    pca = FALSE,
                                    check_duplicates = FALSE
      )
      cost[i] <- tSNE.out[[i]]$costs[which.min(tSNE.out[[i]]$costs)]
    }
  }
  else {
    tSNE.out <- future.apply::future_lapply(seq_len(n.samples), 
                                            function(i) {
      cat("tSNE sample =", i, "of", n.samples, ".\n")
      Rtsne::Rtsne(Z,
                   dims = 2, 
                   perplexity = perp,
                   max_iter = n.iter,
                   verbose = TRUE,
                   pca = FALSE,
                   check_duplicates = FALSE
      )
    }, future.seed = TRUE)
    cost <- unlist(lapply(tSNE.out,
                          function(x) x$costs[which.min(x$costs)]))
  }

  # Getting optimal embedding.
  S.tsne <- tSNE.out[[which.min(cost)]]
  rownames(S.tsne$Y) <- sample.names
  return(S.tsne)
}
