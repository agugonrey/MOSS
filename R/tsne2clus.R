#' t-Stochastic Neighbor Embedding to Clusters
#'
#' Finds clusters on a 2 dimensional map using
#' Density-based spatial clustering of applications with noise
#' (DBSCAN; Esther et al. 1996).
#'
#' The function takes the outcome of pca2tsne (or a list containing any
#'  two-columns matrix) and finds clusters via DBSCAN.
#'  It extends code from the MEREDITH (Taskesen et al. 2016) and
#'  clValid (Datta & Datta, 2018) R packages to tune DBSCAN parameters
#'  with Silhouette or
#'  Biological Homogeneity indexes.
#' @param S.tsne Outcome of function "pca2tsne"
#' @param ann Subjects' annotation data.
#' An incidence matrix assigning subjects to classes of biological
#'  relevance.
#'  Meant to tune cluster assignation via Biological Homogeneity Index (BHI).
#'   If ann=NULL, the number of clusters is tuned with the
#'    Silhouette index instead of BHI. Defaults to NULL.
#' @param labels Character vector with labels describing subjects.
#' Meant to assign aesthetics to the visual display of clusters.
#' @param aest Data frame containing points shape and color.
#' Defaults to NULL.
#' @param eps_range Vector containing the minimum and
#' maximum eps values to be explored. Defaults to c(0, 4).
#' @param eps_res How many eps values should be explored between the 
#' specified range?
#' @param min.clus.size Minimum size for a cluster to appear in the visual
#' display. Defaults to 10
#' @param group.names The title for the legend's key if 'aest' is specified.
#' @param xlab Name of the 'xlab'. Defaults to "x: tSNE(X)"
#' @param ylab Name of the 'ylab'. Defaults to "y: tSNE(X)"
#' @param clus Should we do clustering? Defaults to TRUE. If false, only
#' point aesthetics are applied.
#' @return \itemize{A list with the results of the DBSCAN
#' clustering and (if argument 'plot'=TRUE) the corresponding
#' graphical displays.
#' \item dbscan.res: a list with the results of the (sparse) SVD,
#' containing:
#'  \itemize{
#'    \item cluster: Cluster partition.
#'    \item eps: Optimal eps according to the Silhouette or Biological
#'    Homogeneity indexes criteria.
#'    \item SIL: Maximum peak in the trajectory of the Silhouette index.
#'    \item BHI: Maximum peak in the trajectory of the Biological
#'    Homogeneity index.
#'  }
#' \item clusters.plot: A ggplot object with the clusters' graphical display.
#' }
#' @references \itemize{
#'    \item Ester, Martin, Martin Ester, Hans-Peter Kriegel,
#'    Jorg Sander, and Xiaowei Xu. 1996.
#'    "A Density-Based Algorithm for Discovering Clusters in
#'     Large Spatial Databases with Noise," 226_231.
#'    \item Hahsler, Michael, and Matthew Piekenbrock. 2017.
#'    "Dbscan: Density Based Clustering of Applications with Noise
#'    (DBSCAN) and Related Algorithms."
#'    https://cran.r-project.org/package=dbscan.
#'    \item Datta, Susmita, and Somnath Datta. 2006. Methods for
#'     Evaluating Clustering Algorithms for Gene Expression Data
#'     Using a Reference Set of Functional Classes.
#'     BMC Bioinformatics 7 (1). BioMed Central:397.
#'    \item Taskesen, Erdogan, Sjoerd M. H. Huisman, Ahmed Mahfouz,
#'     Jesse H. Krijthe, Jeroen de Ridder, Anja van de Stolpe,
#'     Erik van den Akker, Wim Verheagh, and Marcel J. T. Reinders. 2016.
#'      Pan-Cancer Subtyping in a 2D-Map Shows Substructures
#'      That Are Driven by Specific Combinations of Molecular
#'       Characteristics. Scientific Reports 6 (1):24949.
#'  }
#'
#' @export
#' @examples
#' \donttest{
#' library(MOSS)
#' library(viridis)
#' library(cluster)
#' library(annotate)
#'
#' # Using the 'iris' data tow show cluster definition via BHI criterion.
#' set.seed(42)
#' data(iris)
#' # Scaling columns.
#' X <- scale(iris[, -5])
#' # Calling pca2tsne to map the three variables onto a 2-D map.
#' Z <- pca2tsne(X, perp = 30, n.samples = 1, n.iter = 1000)
#' # Using 'species' as previous knoledge to identify clusters.
#' ann <- model.matrix(~ -1 + iris[, 5])
#' # Getting clusters.
#' tsne2clus(Z,
#'   ann = ann,
#'   labels = iris[, 5],
#'   aest = aest.f(iris[, 5]),
#'   group.names = "Species",
#'   eps_range = c(0, 3)
#' )
#'
#' # Example of usage within moss.
#' set.seed(43)
#' sim_blocks <- simulate_data()$sim_blocks
#' out <- moss(sim_blocks[-4],
#'   tSNE = TRUE,
#'   cluster = list(eps_range = c(0, 4), eps_res = 100, min_clus_size = 1),
#'   plot = TRUE
#' )
#' out$clus_plot
#' out$clusters_vs_PCs
#' }
tsne2clus <- function(S.tsne, ann = NULL, labels,
                      aest = NULL, eps_res = 100,
                      eps_range = c(0, 4),
                      min.clus.size = 10, group.names = "Groups",
                      xlab = "x: tSNE(X)", ylab = "y: tSNE(X)", clus = TRUE) {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' needs to be installed for graphical displays.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needs to be installed for graphical displays.")
  }

  # Checking if the right packages are present for dbscan.
  if (clus == TRUE) {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' needs to be installed for clustering.")
    }
    if (is.null(ann) == FALSE) {
      if (!requireNamespace("annotate", quietly = TRUE)) {
        stop("Package 'annotate' needs to be installed 
             for tuning clusters based on BHI.")
      }
      if (!requireNamespace("clValid", quietly = TRUE)) {
        stop("Package 'clValid' needs to be installed 
             for tuning clusters based on BHI.")
      }
    }
  }

  S <- S.tsne$Y

  if (is.null(aest)) aest <- aest.f(labels)

  sample.names <- rownames(S)
  if (is.null(sample.names)) sample.names <- as.character(seq_len(nrow(S)))
  n <- length(sample.names)

  # Creating a data.frame to plot the embedding.
  emb.plot <- data.frame(
    ID = sample.names, x = S[, 1], y = S[, 2],
    label = labels
  )

  emb.plot$col <- aest[emb.plot$label, 2]
  emb.plot$label <- factor(emb.plot$label,
    levels = rownames(aest),
    ordered = TRUE
  )
  emb.plot$shape <- as.numeric(aest[emb.plot$label, 1])

  if (clus) {
    if (length(eps_range) == 1) {
      suppressWarnings(S.clus.res <-
        list(
          cluster = dbscan::dbscan(S,
            eps = eps_range,
            MinPts = 5
          )$cluster,
          eps = eps_range
        ))
    } else {
      if (!is.null(ann)) {
        cat("Finding cluster partition maximizing BHI.\n")
        S.fc <- matrix(F, n, ncol(ann), dimnames = list(
          sample.names,
          colnames(ann)
        ))
        tmp <- intersect(sample.names, rownames(ann))
        S.fc[tmp, ] <- ann[tmp, ]
        S.fc <- S.fc[, colSums(S.fc) > 0]
        suppressWarnings(S.clus.res <- dbscan_BHI(S.tsne$Y,
          annot = S.fc,
          eps_res = eps_res,
          eps_range = eps_range
        ))
      }
      else {
        suppressWarnings(S.clus.res <- dbscan_SH(S,
          eps_res = eps_res,
          eps_range = eps_range
        ))
      }
    }

    if (length(unique(S.clus.res$cluster)) > 2) {
      # Creating a data.frame to plot the embedding.
      emb.plot$cluster <- S.clus.res$cluster

      # Removing out of cluster points.
      emb.plot <- emb.plot[emb.plot$cluster != 0, ]

      # Annotating only "big" clusters.
      tmp <- emb.plot$cluster %in%
        names(table(emb.plot$cluster)[table(emb.plot$cluster) >
          min.clus.size])

      emb.plot <- emb.plot[tmp, ]

      emb.plot$cluster <- factor(emb.plot$cluster,
        levels = unique(emb.plot$cluster),
        labels = seq_len(length(unique(emb.plot$cluster)))
      )
      if (length(unique(emb.plot$cluster)) > 1) {
        # Creating a data.frame for cluster labels.
        cluster.label <- as.data.frame(apply(emb.plot[, c("x", "y")],
          2,
          function(x, y) {
            tapply(
              x,
              y,
              mean
            )
          },
          y = emb.plot$cluster
        ))
        cluster.label$label <- as.integer(rownames(cluster.label))
      }
      else {
        emb.plot$cluster <- "1"
        cluster.label <- as.data.frame(rbind(colMeans(emb.plot[, c("x", "y")])))
        cluster.label$label <- "1"
      }
    }
    else {
      emb.plot$cluster <- "1"
      cluster.label <- as.data.frame(rbind(colMeans(emb.plot[, c("x", "y")])))
      cluster.label$label <- "1"
    }
  }

  ylim <- c(
    range(emb.plot$y)[1] + 0.5 * range(emb.plot$y)[1],
    range(emb.plot$y)[2]
  )

  # Plot clusters
  if (clus) {
    emb.plot$cluster <- factor(emb.plot$cluster)
    q <- ggplot2::ggplot(emb.plot, ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::stat_density2d(
        data = emb.plot,
        ggplot2::aes_string(fill = "cluster"),
        alpha = 0.15, geom = "polygon",
        linetype = 0, n = 100,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = emb.plot,
        ggplot2::aes_string(
          x = "x",
          y = "y",
          col = "label",
          shape = "label"
        ),
        size = 1.5, inherit.aes = TRUE
      ) +
      ggplot2::geom_text(
        data = cluster.label,
        ggplot2::aes_string(
          x = "x",
          y = "y",
          label = "label"
        ),
        col = "darkorange", alpha = 0.9, size = 3,
        fontface = "bold", show.legend = FALSE,
        inherit.aes = TRUE
      ) +
      ggthemes::geom_rangeframe() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 8),
        legend.position = c(0, 0),
        legend.justification = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(xlab) +
      ggplot2::scale_y_continuous(ylab, limits = ylim) +
      ggplot2::scale_colour_manual(
        name = group.names,
        labels = rownames(aest),
        values = aest[, 2],
        guide = ggplot2::guide_legend(
          title = group.names,
          nrow = 4
        )
      ) +
      ggplot2::scale_shape_manual(
        name = group.names,
        labels = rownames(aest),
        values = as.integer(aest[, 1]),
        guide = ggplot2::guide_legend(
          title = group.names,
          nrow = 4
        )
      )
    return(list(clusters.plot = q, dbscan.res = S.clus.res))
  }
  else {
    q <- ggplot2::ggplot(emb.plot, ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_point(
        data = emb.plot, ggplot2::aes_string(
          x = "x", y = "y",
          col = "label",
          shape = "label"
        ),
        size = 1.5, inherit.aes = TRUE
      ) +
      ggthemes::geom_rangeframe() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 8),
        legend.position = c(0, 0),
        legend.justification = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(xlab) +
      ggplot2::scale_y_continuous(ylab, limits = ylim) +
      ggplot2::scale_colour_manual(
        name = group.names,
        labels = rownames(aest),
        values = aest[, 2],
        guide = ggplot2::guide_legend(
          title = group.names,
          nrow = 4
        )
      ) +
      ggplot2::scale_shape_manual(
        name = group.names,
        labels = rownames(aest),
        values = as.integer(aest[, 1]),
        guide = ggplot2::guide_legend(
          title = group.names,
          nrow = 4
        )
      )
    return(list(clusters.plot = q))
  }
}
