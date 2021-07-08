#' Returns signatures of features by groups of subjects
#' 
#' This function is meant to used after moss_select.
#' Its main purpose is to visualize how each selected feature (
#' non-zero loading feature) contributes to each group of subjects by 
#' latent dimension.
#' 
#' @param data.blocks A list of omic blocks as provided to moss.
#' @param moss_select.out The output of moss_select.
#' @param clus_lab A vector of same length than number of subjects 
#' with labels used to visualize clusters. Defaults to NULL.
#' @param th Show the th% features of largest loadings (in absolute value). 
#' Default to th=1 (all the features). Numeric.
#' @param only.candidates Should we plot only candidate features? Logical.
#' @param plot Should the results be plotted? Logical. 
#' Defaults to FALSE
#' @param feature.labels List with with features names for each omic.
#' Defaults to NULL.
#' @return Returns a list with 'signatures', and if plot=TRUE,
#' a ggplot object named 'sig_plot'. The element 'signatures'
#' is a data frame with columns corresponding to 'Cluster' (groups
#' of subjects), 'Omic', 'Dim' (PC index or latent dimension),
#' 'Feature_name', 'Feature_pos' (column index of the selected
#' feature within the corresponding omic), 'Loadings' (non-zero 
#' loadings from moss), 'Means', 'L1' and 'L2' (mean +/- standard
#' error of the selected feature values within an omic).
#' @export
#' @examples
#'\donttest{
#' library("MOSS")
#' # Extracting simulated omic blocks.
#' sim_data <- simulate_data()
#' sim_blocks <- sim_data$sim_blocks
#'
#' # Extracting subjects and features labels.
#' lab.sub <- sim_data$labels$lab.sub
#'
#' out <- moss(sim_blocks[-4],
#'   method = "pca",
#'   nu.v = 10,
#'   exact.dg = TRUE,
#'   plot = TRUE,
#'   alpha.v = 0.5
#' )
#' out2 <- moss_select(data.blocks = sim_blocks[-4],
#'                     SVD = out$sparse,
#'                     plot = TRUE)
#'
#' # Display signature plots.
#' out3 <- moss_signatures(data.blocks = sim_blocks[-4],
#'                         clus_lab=lab.sub,
#'                         moss_select.out = out2,
#'                         plot = TRUE)
#' out3$sig_plot 
#' }                               
moss_signatures <- function(data.blocks, 
                            moss_select.out, 
                            clus_lab =NULL,
                            plot=FALSE, 
                            feature.labels=NULL,
                            th=1,
                            only.candidates=FALSE) {

  M <- length(data.blocks)
  
  # Naming data blocks if necessary.
  if (is.null(names(data.blocks))) {
    names(data.blocks) <- paste(
      "Block",
      seq_len(M)
    )
  } 
  else {
    tmp <- names(data.blocks) == ""
    names(data.blocks)[tmp] <- paste("Block", seq_len(M))[tmp]
  }
  
  # Naming features if necessary.
  if (is.null(feature.labels)) {
    feature.labels <- lapply(data.blocks, function(x) {
      if (is.null(colnames(x))) {
        paste(
          "Feature",
          seq_len(ncol(x)))
      } 
      else colnames(x)
    })  
  }
  names(feature.labels) <- names(data.blocks)

  tmp <- which(!(grepl(names(moss_select.out), pattern = "Subjects")| 
      grepl(names(moss_select.out), pattern = "features_contr")))
  omic_names <- names(moss_select.out)[tmp]
  feature.labels <- feature.labels[omic_names]
  K <- names(moss_select.out[[tmp[1]]])
  signatures <- NULL
  
  if (is.null(clus_lab)) clus_lab <- rep(1,nrow(data.blocks[[1]]))
  for (cl in unique(clus_lab)) {
    i <- which(clus_lab == cl)
    for (r in tmp) {
      for (k in K) {
        tmp2 <- as.numeric(names(moss_select.out[[r]][[k]]))
        if (length(tmp2) > 0) {
          sig <- data.frame("Cluster" = cl,
                            "Omic" = omic_names[r],
                            "Dim" = k,
                          "Feature_name" = feature.labels[[r]][tmp2],
                            "Feature_pos" = tmp2,
                            "Loadings" = moss_select.out[[r]][[k]])
          sig$Omic <- factor(sig$Omic,
                             levels=omic_names,ordered = TRUE)
          aux <-  as.data.frame(do.call("rbind",
                                lapply(tmp2, 
                                  function(j) {
                                    x <- data.blocks[[r]][i, j]
                                    mu <- mean(x, na.rm = TRUE)
                                    s <- stats::sd(x, na.rm = TRUE) /
                                    sqrt(length(mu))
                                    return(c(mu - s, mu, mu +s))
                                               }
                                        )))
          names(aux) <- c("L1","Means","L2")
          sig <- cbind(sig, aux)
          aux <- data.frame("Candidate"=apply(sign(sig[,c("L1","L2")]),
                                  1,
                                  prod) > 0)
          sig <- cbind(sig, aux)
          signatures <- rbind(signatures, sig)
        } 
      }
    }
  }
  out <- NULL
  out$signatures <- signatures
  
  if (only.candidates) signatures <- 
    signatures[signatures$Candidate == TRUE,]
  
  signatures <- signatures[signatures$Loadings ^ 2 >= 
    stats::quantile(signatures$Loadings ^ 2,probs = 1 - th), ]
  
  if (plot) {
    g <- vector("list",length(unique(clus_lab)))
    names(g) <- unique(clus_lab)
    for (cl in unique(clus_lab)) {
      g[[as.character(cl)]] <- 
        ggplot2::ggplot(signatures[signatures$Cluster == cl,],
                                 ggplot2::aes_string(
                                   x = "Means",
                                   y = "Feature_name",col="Omic")) +
        ggplot2::geom_point()+
        ggplot2::geom_errorbarh(ggplot2::aes_string(xmax = "L2", 
                                           xmin = "L1",
                                           height = .2))+
        ggplot2::facet_grid(stats::as.formula("Omic ~ Dim"),
                            scales="free")+
ggplot2::scale_color_manual(values = viridis::viridis(length(tmp)))+
        ggplot2::scale_x_continuous("Mean (SE)") +
        ggplot2::scale_y_discrete("Feature name") +
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "none")+
        ggplot2::geom_vline(xintercept = 0,
                            linetype="dotted",
                            color="#00000080")
    }
    out$sig_plot <- g
  }
  return(out)
}