#' Simple simulation of regulatory modules across three omic blocks and one classification variable.
#'
#' This a simple simulation to use in MOSS' examples. The specifics of the simulation are shown in the "Examples" section.
#' @param moss_seed The seed for random number generator. Numeric. Defaults to 42.
#' @return
#' A list of two elements 'sim_blocks' and 'labels'. 
#' First element 'sim_blocks' is a list of three numeric matrices, and one character matrix. Second element 'labels' has two character vectors. 
#' The first element 'lab.sub' identifies the groups of 'signal' subjects.
#' The second element 'lab.feat' identifies the groups 'signal' features from background 'noise'.
#' @export
#' @examples
#' sim_data <- simulate_data()
#'
#' #Extracting simulated omic blocks.
#' sim_blocks <- sim_data$sim_blocks
#' 
#' #Extracting subjects and features labels.
#' lab.sub <- sim_data$labels$lab.sub
#' lab.feat <- sim_data$labels$lab.feat
#'
#' #Check dimensions and objects class.
#' lapply(sim_blocks,dim)
#' lapply(sim_blocks,function(x) class(x[,1]))
#'
#' #Showing how the data was generated.
#' set.seed(42)
#' O1 <- matrix(data=0, nrow=5e2,ncol=1e3)
#' O2 <- O1
#' O1[1:20, 1:150] <- 1
#' O1 <- O1 + rnorm(5e5, mean=0,sd=0.5)
#' O2[71:90, 71:200] <- 1
#' O2 <- O2 + rnorm(5e5, mean=0,sd=0.5)
#' #Simulating a continous response blocks.
#' O3 <- 3*O1 - 5*O2 + rnorm(5e5, mean=0,sd=0.5)
#' 
#' #Creating a vector labeling clusters of subjects.
#' aux <- rep("Background",500)
#' aux[1:20] <- "Group 1"
#' aux[71:90] <- "Group 2"
#' all.equal(aux,lab.sub)
#' 
#' #Generating a classification response.
#' O4 <- as.matrix(aux)
#' 
#' #Storing all blocks within a list.
#' all.equal(sim_blocks,list("Block 1" = O1,
#'    "Block 2" = O2,
#'    "Block 3" = O3,
#'    "Block 4" = O4))
#'
#' #Creating a vector labeling signal and background features.
#' aux <- rep("Background features",3000)
#' aux[c(1:150,1072:1200,2001:2200)] <- "Signal features"
#' all.equal(aux,lab.feat)
simulate_data <- function(moss_seed = 42) {
  #Simple simulation of three omic blocks and a single categorical variable.
  set.seed(moss_seed)
  O1 <- matrix(data=0, nrow=5e2,ncol=1e3)
  O2 <- O1
  O1[1:20, 1:150] <- 1
  O1 <- O1 + stats::rnorm(5e5, mean=0,sd=0.5)
  
  O2[71:90, 71:200] <- 1
  O2 <- O2 + stats::rnorm(5e5, mean=0,sd=0.5)
  
  #Simulating a continous response blocks.
  O3 <- 3*O1 - 5*O2 + stats::rnorm(5e5, mean=0,sd=0.5)
  
  #Creating a vector labeling clusters of subjects.
  lab.sub <- rep("Background",500)
  lab.sub[1:20] <- "Group 1"
  lab.sub[71:90] <- "Group 2"
  
  #Generating a classification response.
  O4 <- as.matrix(lab.sub)
  
  #Storing all blocks within a list.
  sim_blocks <- list("Block 1" = O1,
                     "Block 2" = O2,
                     "Block 3" = O3,
                     "Block 4" = O4)
  
  #Creating a vector labeling signal and background features.
  lab.feat <- rep("Background features",3000)
  lab.feat[c(1:150,1072:1200,2001:2200)] <- "Signal features"
  
  #Storing simulated data in a global variable.
  
  sim_data <- list("sim_blocks"=sim_blocks,
                    "labels"=list("lab.sub"=lab.sub,
                                  "lab.feat"=lab.feat)) 
  return(sim_data)
}
