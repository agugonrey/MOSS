#' Extracts (and merges) chunks of characters.
#'
#' @param x A character vector.
#' @param i Index specifying which chunks of characters 
#' will be extracted (and merged).
#' @param sep Chunks separator character. Defaults to "-".
#' @param collapse New chunks separator character.
#' Default to 'sep'.
#' @returns A character vector with the extracted (and merged) chunks of
#' characters.
#' @export
#' @examples
#' x <- "this is one chunk of characters & this is another one"
#' metdat(x, 1, " & ")
#' metdat(x, 2, " & ")
#' metdat(x, c(1, 2), " & ")
#' metdat(x, c(1, 2), " & ", " and ")
metdat <- function(x, i, sep = "-", collapse = sep) {
  unlist(lapply(strsplit(x, sep, T), function(x) {
    paste(x[i],
      collapse = collapse
    )
  }))
}
