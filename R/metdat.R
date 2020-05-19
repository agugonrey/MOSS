metdat <- function(x,i,sep = "-") unlist(lapply(strsplit(x,sep,T), function(x) paste(x[i], collapse = "-")))
