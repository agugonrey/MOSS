# A function to calculate the Biological Homogeneity Index (from clusteval package)
BHI <- function (statClust, annotation, names = NULL, category = "all", dropEvidence = NULL) {
  if (is.matrix(annotation)) {
    bhi <- numeric(length(unique(statClust)))
    names(bhi) <- unique(statClust)
    for (k in unique(statClust)) {
      Ck.bhi <- 0
      Ck.idx <- which(statClust == k)
      if (length(Ck.idx) < 2) 
        next
      for (i in Ck.idx) {
        B <- which(annotation[i, ] == TRUE)
        if (length(B) == 0) 
          next
        annot <- annotation[Ck.idx[Ck.idx != i], B]
        if (length(B) == 1) 
          Ck.bhi <- Ck.bhi + sum(annot)
        else if (length(B) > 1) 
          if (is.null(dim(annot))) Ck.bhi <- Ck.bhi + sum(sum(annot) > 0)
          else Ck.bhi <- Ck.bhi + sum(rowSums(annot) > 0)
      }
      nk <- sum(rowSums(annotation[Ck.idx, ]) > 0)
      if (nk > 1) 
        bhi[k] <- Ck.bhi/(nk * (nk - 1))
    }
    return(mean(bhi, na.rm = TRUE))
  }
  goTerms <- annotate::getGO(names, annotation)
  if (!is.null(dropEvidence)) 
    goTerms <- lapply(goTerms, annotate::dropECode, dropEvidence)
  bhi <- tapply(goTerms, statClust, function(x) clValid::matchGO(x, category))
  bhi.mean <- mean(bhi, na.rm = TRUE)
  bhi.se <- stats::sd(bhi, na.rm = TRUE) / sqrt(length(bhi))
  return(bhi.mean)
}
