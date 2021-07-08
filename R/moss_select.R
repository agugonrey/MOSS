#' Returns features and subject selected by latent dimension.
#' 
#' This function is meant to used after moss.
#' Its main purpose is to extract the features and subjects by 
#' latent dimension.
#' The selection depends on loadings at each dimension being 
#' different from zero.
#' 
#' @param data.blocks a list of omic blocks as provided to moss.
#' @param SVD a list with SVD results. The function is meant to work
#'  with the results from sparse SVD. 
#'  However, 'dense' solutions are also accepted.
#' @param resp.block Which omic block was used as response in moss? 
#' Integer. Defaults to NULL.
#' @param K How many dimensions should be displayed? Vector.
#' Defaults to the 1 : ncol(SVD$v).
#' @param plot Should the results be plotted? Logical. 
#' Defaults to FALSE
#' @return Returns a list containing the position, name, and loadings
#' of selected features and subjects by latent dimension.
#' if 'plot=TRUE', a scatterplot is displayed, where the x-axis
#' represents the latent dimensions, the y-axis the total number 
#' of features selected in log scale, and each point is a pie chart
#' showing the relative contribution of each omic to the number of
#' features selected. The radio of the pie chart represents the 
#' coefficient of variation among squared loadings 
#' (mean squared loadings divided by their standard deviation).
#' @export
moss_select <- function(data.blocks, 
                       SVD, 
                       resp.block = NULL,
                       K=NULL,
                       plot=FALSE) {
  
  SVD$u <- as.matrix(SVD$u)
  SVD$v <- as.matrix(SVD$v)
  
  if (!is.null(K)) {
    SVD$u <- SVD$u[, K]
    SVD$v <- SVD$v[, K]
  }
  else K <- seq_len(ncol(SVD$u))
  
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

  if (is.null(resp.block)) {
    p <- lapply(data.blocks, ncol)
    j <- cumsum(p)
    j <- lapply(X = seq_len(length(j)),
                FUN = function(i) {
                  s <- seq(ifelse(length(j[i - 1]) == 0,
                                  0,
                                  j[i - 1]) + 1,
                           j[i])
                  names(s) <- seq_along(s)
                  s
                })
    omic_indexes <- unlist(lapply(j,names))
    omic_blocks <- seq_len(length(data.blocks))
    feat_sel <- vector("list", length = length(omic_blocks))
    names(feat_sel) <- names(data.blocks)[omic_blocks]
    for (k in seq_along(omic_blocks)) {
      feat_sel[[k]] <- lapply(seq_len(ncol(SVD$v)), 
                              function(i) {
                                n <- which(SVD$v[j[[k]],i] != 0)
                                v <- SVD$v[j[[k]],i][n]
                                names(v) <- omic_indexes[n]
                                return(v)
                                })
      
      names(feat_sel[[k]]) <- paste0("Dim",K)
    }
    feat_sel[[length(data.blocks)+1]] <- 
      lapply(seq_len(ncol(SVD$u)),
             function(i) {
               n <- which(SVD$u[,i] != 0)
               v <- SVD$u[,i][n]
               names(v) <- n
               return(v)
             })
    names(feat_sel)[length(data.blocks)+1] <- "Subjects"
  }
  
  else {
    p <- lapply(data.blocks[-resp.block], ncol)
    j <- cumsum(p)
    j <- lapply(X = seq_len(length(j)),
           FUN = function(i) {
             s <- seq(ifelse(length(j[i - 1]) == 0,
                                        0,
                                        j[i - 1]) + 1,
                                 j[i])
             names(s) <- seq_along(s)
             s
             })
    omic_indexes <- unlist(lapply(j,names))
    
    omic_blocks <- seq_len(length(data.blocks))[-resp.block]
    feat_sel <- vector("list", length = length(omic_blocks))
    names(feat_sel) <- names(data.blocks)[omic_blocks]
    for (k in seq_along(omic_blocks)) {
      feat_sel[[k]] <- lapply(seq_len(ncol(SVD$u)), 
                              function(i) {
                                n <- which(SVD$u[j[[k]],i] != 0)
                                v <- SVD$u[j[[k]],i][n]
                                names(v) <- omic_indexes[n]
                                return(v)
                              })
      names(feat_sel[[k]]) <- paste0("Dim",seq_len(ncol(SVD$u)))
    }
    feat_sel[[length(data.blocks)]] <- 
      lapply(seq_len(ncol(SVD$v)),
             function(i) {
               n <- which(SVD$v[,i] != 0)
               v <- SVD$v[,i][n]
               names(v) <- n
               return(v)
             })
    names(feat_sel[[length(data.blocks)]]) <- 
      paste0("Dim",K)
    
    names(feat_sel)[length(data.blocks)] <- 
      names(data.blocks)[resp.block]
    feat_sel <- c(feat_sel[length(data.blocks)],
                  feat_sel[-length(data.blocks)])
    
  }
  
  if (plot == TRUE & length(data.blocks) > 1) {
  
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needs to be installed for graphical 
           displays.")
      }
     
      if (!requireNamespace("scatterpie", quietly = TRUE)) {
    stop("Package 'scatterpie' needs to be installed for graphical 
           displays.")}
        
      if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' needs to be installed for graphical 
           displays.")
      }
    
    aux <- feat_sel[names(feat_sel) != "Subjects"]
    d <- list("x" = seq_len(ncol(SVD$u)))

    for (k in seq_along(aux)) {
      d[[k + 1]] <- rep(0,ncol(SVD$u))
      for (i in seq_len(ncol(SVD$u))) {
        d[[k + 1]][i] <- length(aux[[k]][[i]])
      }
    }
    names(d)[names(d) != "x"] <- names(aux)
    d$y <- seq_len(ncol(SVD$u))
    d$radius <- seq_len(ncol(SVD$u))
    for (i in seq_len(ncol(SVD$u))) {
      v <- NULL
      for (k in seq_along(aux)) {
        v <- c(v, aux[[k]][[i]])
      }
      d$y[i] <- mean(v) / stats::sd(v)
      d$radius[i] <- length(v)
    }
   # d$radius <- (d$radius - min(d$radius)) / diff(range(d$radius))
    d$radius <- log(d$radius,2)
    d$y <- 0.15 * (d$y - min(d$y)) / diff(range(d$y)) + 0.05
    d <- as.data.frame(do.call("cbind",d))
    r <- tapply(d$radius,d$x,max)
    a <- which.min(r); names(a) <- NULL
    b <- max(r) + d$y[a]; names(b) <- NULL
    g <- ggplot2::ggplot(d,
      ggplot2::aes_string(
        x = "x",
        y = "radius")) +
      ggplot2::geom_point()+
      scatterpie::geom_scatterpie(ggplot2::aes_string(x="x", 
                                                      y="radius", 
                                                      group="x", 
                                                      r="y"),
                                  data=d, pie_scale =T ,
                                  cols=names(aux), 
                                  color=NA, alpha=.8) +
      scatterpie::geom_scatterpie_legend(round(d$y,2),
                                         x = ncol(SVD$u)-0.5,y=b)+
  ggplot2::scale_fill_manual(values = viridis::viridis(length(aux)))+
      ggplot2::scale_x_continuous("PC index",
                                  breaks = seq_len(ncol(SVD$u)), 
                                  labels = K) +
      ggplot2::scale_y_continuous("log(Number of features)") +
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.position = "top")+
      ggplot2::guides("fill"=ggplot2::guide_legend(title = "Omic blocks",
                                          title.position = "top",
                                          title.hjust = 0.5))
     
    feat_sel$features_contr <- g 
  }
      
  return(feat_sel)
}