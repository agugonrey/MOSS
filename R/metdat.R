#' Extracts (and merges) characters chunks.
#' 
#' @param x A character vector.
#' @param i Index specifying which one of the character chunks separated by 'sep' is to be picked.
#' @param sep Character separating chunks of characters. Defaults to "-".
#' @param collapse New character vector separating chunks of characters. Default to 'sep'.
#' @returns A character vector with the extracted (and merged) chunks of characters from input vector.
#' @export
#' @examples
#' x <- "this is one chunk of characters & this is another" 
#' metdat(x, 1, " & ")
#' metdat(x, 2, " & ")
#' metdat(x, c(1, 2), " & ")
#' metdat(x, c(1, 2), " & ", " and ")
metdat <- function(x,i,sep = "-", collapse = sep) unlist(lapply(strsplit(x,sep,T), function(x) paste(x[i], collapse = collapse)))
