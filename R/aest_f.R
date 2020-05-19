aest.f <- function(x,n.cat=2,option="D") {
  if(!requireNamespace("viridis",quietly = TRUE)) stop("Package 'viridis' needs to be installed.")
  if (is.numeric(x) == T) x <- ggplot2::cut_interval(x, n.cat)
  if (is.factor(x) == T) x <- as.character(x)
  d <- data.frame(pch=16, col=viridis::viridis(length(unique(x)),option = option,alpha=0.5),stringsAsFactors = F)
  names <- unique(x)
  tmp <- is.na(names)
  d <- d[!tmp,]
  rownames(d) <- names[!tmp]
  return(d)
}
