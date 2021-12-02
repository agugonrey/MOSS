#' Useful Venn diagrams to study the overlap between samples row names.
#'
#' @param L List of elements which overlap we wish to check (e.g., row names 
#' by omic blocks).
#' @param a Elements of the list we want to focus on (e.g., a subset of omic 
#' blocks). Numerical.
#' @param lty Line width of circles circumferences.
#' @param fill Color for each circle. Character vector. Defaults to NULL.
#' @param element_names Names of each category. Character vector.
#' Defaults to NULL
#' @return A data.frame with labels as rownames and two
#' @export
moss_venn <- function(L,a,lty="blank",fill=NULL,element_names=NULL) {
  areas <- function(L,i,areas=T) {
    ints <- L[[i[1]]]
    if (length(i) > 1) {
      for(l in i[-1]) {
        ints <- intersect(ints, L[[l]])
      }
    }
    if (areas) return(length(ints))
    else return(ints)
  }
  if (is.null(element_names)) element_names <- seq_along(L)
  if (is.null(names(L))) names(L) <- element_names
  if (is.null(fill)) fill <- viridis::viridis(length(a),alpha = 0.7)
  category <- names(L)
  grid::grid.newpage()
  if (length(a) == 1) {
    out <- VennDiagram::draw.single.venn(areas(L,a), category=category,lty=lty,fill=fill)
  }
  if (length(a) == 2) {
    out <- VennDiagram::draw.pairwise.venn(area1 = areas(L, a[1]), areas(L,a[2]), areas(L,a[1:2]), category=category,lty=lty,fill=fill)
  }
  if (length(a) == 3) {
    out <- VennDiagram::draw.triple.venn(areas(L, a[1]), areas(L,a[2]), areas(L,a[3]), areas(L,a[1:2]), 
                                         areas(L,a[2:3]), areas(L,a[c(1, 3)]), areas(L,a), category=category,lty=lty,fill=fill)
  }
  if (length(a) == 4) {
    out <- VennDiagram::draw.quad.venn(areas(L,a[1]), areas(L,a[2]), areas(L,a[3]), areas(L,a[4]), 
                                       areas(L,a[1:2]), areas(L,a[c(1, 3)]), areas(L,a[c(1, 4)]), areas(L,a[2:3]), 
                                       areas(L,a[c(2, 4)]), areas(L,a[3:4]), areas(L,a[1:3]), areas(L,a[c(1, 2, 
                                                                                                          4)]), areas(L,a[c(1, 3, 4)]), areas(L,a[2:4]), areas(L,a), category=category,lty=lty,fill=fill)
  }
  if (length(a) == 5) {
    out <- VennDiagram::draw.quintuple.venn(areas(L,a[1]), areas(L,a[2]), areas(L,a[3]), areas(L,a[4]), areas(L,a[5]),
                                            areas(L,a[1:2]), areas(L,a[c(1, 3)]), areas(L,a[c(1, 4)]),areas(L,a[c(1, 5)]),
                                            areas(L,a[2:3]), areas(L,a[c(2, 4)]), areas(L,a[c(2, 5)]), 
                                            areas(L,a[3:4]),areas(L,a[c(3,5)]), areas(L,a[c(4,5)]), areas(L,a[1:3]), areas(L,a[c(1, 2, 4)]), areas(L,a[c(1, 2, 5)]),
                                            areas(L,a[c(1, 3, 4)]), areas(L,a[c(1, 3, 5)]), areas(L,a[c(1, 4, 5)]),areas(L,a[2:4]), 
                                            areas(L,a[c(2, 3, 5)]), areas(L,a[c(2, 4, 5)]), areas(L,a[c(3, 4, 5)]), areas(L,a[c(1:4)]),
                                            areas(L,a[c(1, 2, 3, 5)]),areas(L,a[c(1, 2, 4, 5)]),areas(L,a[c(1, 3, 4, 5)]),
                                            areas(L,a[c(2, 3, 4, 5)]),areas(L,a[c(1:5)]), areas(L,a), category=category,lty=lty,fill=fill)
  }
  if (!exists("out")) warning("Could not produce plot :( !")
  grid::grid.draw(out)
}