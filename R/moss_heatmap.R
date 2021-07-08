#' Creates a heatmap from the output of MOSS.
#' 
#' @param B An bject of class 'matrix' or 'FBM'.
#' @param SVD List with the results of a sparse or dense SVD.
#' @param right.lab Columns title. Character.
#' @param left.lab Rows title. Character.
#' @param axes.pos What SVD dimensions should be used to plot the heatmap? 
#' If NULL, all the SVD dimensions are use. Defaults to NULL.
#' @param verbose Should we print messages? Logical.
#' Defaults to TRUE.
#' @return Returns a 'ComplexHeatmap' plot representing the cross-product
#'  between left and right Eigenvectors.
#' @export
moss_heatmap <- function(B,SVD,right.lab, left.lab,axes.pos=NULL,verbose=TRUE) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' needs to be installed 
           to generate heatmaps.")
  }  
  else {
    
    if (!is.null(axes.pos)) {
      SVD$u <- SVD$u[,axes.pos]
      SVD$v <- SVD$v[,axes.pos]
    }
    
    if (verbose) message("Creating heatmap.")
    B.sub <- B[which(x = rowMeans(x = SVD$u[,axes.pos] != 0) > 0),
                   which(rowMeans(SVD$v[,axes.pos] != 0) > 0)]
    
    heat_plot <- ComplexHeatmap::Heatmap(matrix = B.sub,
                                             column_title = right.lab,
                                             row_title = left.lab,
                                             show_column_names = FALSE,
                                             show_row_names = FALSE,
                                             show_heatmap_legend = FALSE)
    return(heat_plot)
  }
}