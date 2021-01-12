#' Assign point color and shape aesthetics.
#'
#' This function is called by moss whenever a plot is produced.
#' It simply assigns colors and shape to points based on input labels.
#' @param x Character vector with labels, or a numerical vector to be 
#' discretized in 'n.cat' categories.
#' @param n.cat Number of categories to split vector 'x'. Numeric.
#' Ignored if 'x' is a character vector.
#' @param option Controls color palette.
#' One of the possible 'option' arguments for the 'viridis' function.
#' @return A data.frame with labels as rownames and two
#' columns representing point colors and shape, respectively.
#' @export
aest.f <- function(x, n.cat = 2, option = "D") {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' needs to be installed.")
  }
  if (is.numeric(x) == T) x <- ggplot2::cut_interval(x, n.cat)
  if (is.factor(x) == T) x <- as.character(x)
  d <- data.frame(
    pch = 16,
    col = viridis::viridis(length(unique(x)),
      option = option, alpha = 0.5
    ),
    stringsAsFactors = F
  )
  names <- unique(x)
  tmp <- is.na(names)
  d <- d[!tmp, ]
  rownames(d) <- names[!tmp]
  return(d)
}
