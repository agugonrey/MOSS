#' Multi-omic integration via singular value decomposition.
#'
#' This function concatenates multiple omic, allowing one of them to be a multivariate numeric response 'Y', or a univariate classification response (to allow both unsupervised and supervised omic integration).
#' In general, omic blocks consisting of predictors are concatenated and normalized to create an extended omic matrix 'Z'. Generalized SVD (basically, a SVD on a function of the cross-product between two matrices) allows to integrate omics according to many different multi-variate techniques.
#' All the different applications of multivariate techniques within MOSS return a matrix 'B'. MOSS allows to fit many different multivariate techniques (e.g. B = Z in pca; B = Z'Y, for pls; B = (Z'Z)^-Z'Y, for rrr).
#'
#' Once 'dense' solutions for each technique are found (the result of SVD on a matrix B), the function ssvdEN_sol_path is called to perform sparse SVD (sSVD) on a grid of possible degrees of sparsity.
#' The sSVD is performed using the algorithm of Shen and Huang (2008), extended to include Elastic Net type of regularization. For one latent factor (rank 1 case), the algorithm finds vectors Q, lambda, and T that minimize:
#'     \tabular{rl}{
#'  \tab  ||B - lambda * QT||_F^2 + delta_T(alpha_T||T||_1 + (1 - alpha_T)||T||_F^2) + delta_Q (alpha_Q||Q||_1 + (1 - alpha_Q)||Q||_F^2) \cr
#' }
#'     
#' such that ||Q|| = 1. The right Eigen vector is obtained from T / ||T|| and the corresponding lambda = Q'BT.
#'
#' The penalties delta_T and delta_Q are mapped from specified desired degree of sparsity.
#' Selecting degree of sparsity: The function allows to tune the degree of sparsity using an ad-hoc method based on the one presented in Shen & Huang (2008, see reference) and generalized for tuning two sparsity degrees parameters.
#' This is done by exploring the proportion of explained variance (PEV) at each value of a grid of possible values.
#' Drastic and/or steep changes in the PEV trajectory accross degrees of sparsity are used for automatical selection (see help for the function ssvdEN_sol_path).
#' By impossing the additional assumption of omic blocks being conditionally independent, each multivariate technique can be extended using a 'multi-block' approach, where the contribution of each omic block to the total (co)variance is addressed.
#' When response Y is a character column matrix, with classes or categories by subject, each multivariate technique can be extended to perform linear discriminant analysis.
#'
#' @note
#' \enumerate{
#'   \item The function does not return PEV for EN parameter (alpha-T and/or alpha_Q), the user needs to provide a single value for each.
#'   \item When number of PC index > 1, columns of T might not be orthogonal.
#'   \item Although the user is encouraged to perform data projection and cluster separately, MOSS allows to do this automatically. However, both tasks might require finner tuning than the provided by default, and computations could become cumbersome.
#'   \item Tuning of degrees of sparsity is done heuristically on training set. In our experience, this results in high specificity, but rather low sensitiviy (i.e. too liberal cutoffs, as compared with extensive cross-validation on testing set).
#'   \item When 'method' is an unsupervised technique, 'K.X' is the number of latent factors returned and used in further analysis. When 'method' is a supervised technique, 'K.X' is used to perform a SVD to facilitate the product of large matrices and inverses.
#'   \item If 'K.X' (or 'K.Y') equal 1, no plots are returned.
#'   \item Although the degree of sparsity maps onto number of features/subjects for Lasso, the user needs to be aware that this conceptual correspondence
#'         is lost for full EN (alpha belonging to (0, 1); e.g. the number of features selected with alpha < 1 will be eventually larger than the optimal degree of sparsity).
#'         This allows to rapidly increase the number of non-zero elements when tuning the degrees of sparsity. 
#'         In order to get exact values for the degrees of sparsity at subjects or features levels, the user needs to 
#'         set the value of 'exact.dg' parameter from 'FALSE' (the default) to 'TRUE'.
#' }
#' @param data.blocks List containing omic blocks of class 'matrix' or 'FBM'. In each block, rows represent subjects and columns features. IMPORTANT: omic blocks have to be aligned by rows.
#' @param method Multivariate method. Character. Defaults to 'pca'. Possible options are pca, mbpca, pca-lda, mbpca-lda, pls, mbpls, pls-lda, mbpls-lda, rrr, mbrrr, rrr-lda, mbrrr-lda.
#' @param resp.block What block should be used as response? Integer. Only used when the specified method is supervised.
#' @param K.X Number of principal components for predictors. Integer. Defaults to 5.
#' @param K.Y Number of responses PC index when method is supervised. Defaults to K.X.
#' @param verbose Should we print messages? Logical. Defaults to TRUE.
#' @param ncores Number of cores used for sSVD. Only relevant when at least one omic block is a FBM. Defaults to 1.
#' @param dg.grid.left A grid with increasing integers representing degrees of sparsity for left-eigenvectors. Defaults to NULL.
#' @param dg.grid.right Same but for right eigen vectors. Defaults to NULL.
#' @param alpha.right Elastic Net parameter for right eigenvectors. Numeric between 0 and 1. Defaults to 1.
#' @param alpha.left  Elastic Net parameter for right eigenvectors. Numeric between 0 and 1. Defaults to 1.
#' @param clus Should cluster be obtained? Logical. Defaults to FALSE. If TRUE, DBSCAN is called.
#' @param clus.lab A vector of same length than number of subjects with labels used to visualize clusters. Factor. Defaults to NULL. 
#' When sparsity is imposed on the left eigenvectors, the association between non-zero loadings and labels' groups is shown by a Chi-2 statistics for each pc. When sparsity is not imposed, the association between labels and PC is addressed by a Kruskal-Wallis statistics.
#' @param tSNE Arguments passed to the function pca2tsne as a list. Defaults to NULL. If tSNE=T, defaults parameters are used (perp=50,n.samples=1,n.iter=1e3).
#' @param axes.pos PC index used for tSNE. Defaults to 1 : K.Y. Used only when tSNE is different than NULL.
#' @param scale.arg Should the omic blocks be centered and scaled? Logical. Defaults to TRUE.
#' @param norm.arg Should omic blocks be normalized? Logical. Defaults to TRUE.
#' @param plot Should results be plotted? Logical. Defaults to FALSE.
#' @param approx.arg Should we use standard SVD or random approximations? Defaults to FALSE. If TRUE and at least one block is of class 'matrix', irlba is called. If TRUE & is(O,'FBM')==TRUE, big_randomSVD is called.
#' @param exact.dg Should we compute exact degrees of sparsity? Logical. Defaults to FALSE. Only relevant When alpha.s or alpha.f are in the (0,1) interval and exact.dg = TRUE.
#' @param use_fbm Should we treat omic blocks as Filed Backed Matrix (FBM)? Logical. Defaults to FALSE.

#' @return Returns a list with the results of the sparse generalized SVD. If \emph{plot}=TRUE, a series of plots is generated as well.
#' \itemize{
#' \item \emph{\strong{B:}}  The object of the (sparse) SVD. Depending of the method used, B can be a extended matrix of normalized omic blocks, a variance-covariance matrix, or a matrix of regression coeficients.
#' If at least one of the blocks in 'data.blocks' is of class FBM, is(B,'FBM') is TRUE. Otherwise, is(B,'matrix') is TRUE.
#' \item \emph{\strong{dense:}} A list containing the resuls of the dense SVD.\itemize{
#'    \item \strong{u:} Matrix with left eigenvectors.
#'    \item \strong{v:} Matrix with right eigenvectors.
#'    \item \strong{d:} Matrix with singular values.
#'  }
#'  \item \emph{\strong{sparse:}} A list containing the results of the sparse SVD.\itemize{
#'    \item \strong{u:} Matrix with left eigenvectors.
#'    \item \strong{v:} Matrix with right eigenvectors.
#'    \item \strong{d:} Matrix with singular values.
#'    \item \strong{opt.dg.right:} Selected degrees of sparsity for right eigenvectors. 
#'    \item \strong{opt.dg.left:} Selected degrees of sparsity for left eigenvectors. 
#'  }
#'  \item Graphical displays: Depending on the values in 'plot','tSNE','clus', and 'clus.lab' arguments, the following ggplot objects can be obtained. They contain:\itemize{
#'    \item \strong{scree.plot:} Plots of eigenvalues and their first and second order empirical derivatives along PC indexes. 
#'    \item \strong{tun_dgSpar.plot:} Plots with the PEV trajectory, as well as its first and second empirical derivatives along the degrees of sparsity path.
#'    \item \strong{PC1_2.plot:} Plot of the first two principal components.
#'    \item \strong{tSNE.plot:} Plot with the tSNE mapping onto two dimensions.
#'    \item \strong{clus.obj:} The output of function tsne2clus.
#'    \item \strong{subLabels_vs_cluster:} Plot of the Kruskal-Wallis (or Chi-square) statistics of the association test between PC and pre-established subjects groups.
#'  }
#'  
#' }
#' @references \itemize{
#'    \item Shen, Haipeng, and Jianhua Z. Huang. 2008. Sparse Principal Component Analysis via Regularized Low Rank Matrix approximation. Journal of Multivariate Analysis 99 (6). Academic Press: 1015_34. 
#'    \item Baglama, Jim, Lothar Reichel, and B W Lewis. 2018. Irlba: Fast Truncated Singular Value Decomposition and Principal Components Analysis for Large Dense and Sparse Matrices.
#'  }
#' @export
#' @examples
#' #Example1: sparse PCA of a list of omic blocks.
#' library("MOSS")
#' sim_data <- simulate_data()
#' set.seed(43)
#'
#' #Extracting simulated omic blocks.
#' sim_blocks <- sim_data$sim_blocks
#' 
#' #Extracting subjects and features labels.
#' lab.sub <- sim_data$labels$lab.sub
#' lab.feat <- sim_data$labels$lab.feat
#' out <- moss(sim_blocks[-4],
#'      method="pca",
#'      dg.grid.right = seq(1,200,by=10),
#'      dg.grid.left = seq(1,100,by=2),
#'      alpha.right = 0.5,
#'      alpha.left = 1)
#' 
#' \dontrun{
#' library(ggplot2)
#' library(dbscan)
#' library(ggthemes)
#' library(viridis)
#' library(cluster)
#' library(fpc)
#' 
#' set.seed(43)
#' 
#' #Extracting simulated omic blocks.
#' 
#' #Example2: sparse PCA with t-SNE, clustering, and association with predefined groups of subjects.
#' out <- moss(sim_blocks[-4],
#'      method="pca",
#'      dg.grid.right = seq(1,200,by=10),
#'      dg.grid.left = seq(1,100,by=2),
#'      alpha.right = 0.5,
#'      alpha.left = 1,
#'      tSNE=TRUE,
#'      clus=TRUE,
#'      clus.lab=lab.sub,
#'      plot=TRUE)
#'      
#' #This shows obtained clusters with labels from pre-defined groups of subjects.
#' out$clus.obj
#' 
#' #This shows the statistical overlap between PCs and the pre-defined groups of subjects.
#' out$subLabels_vs_cluster
#' 
#' #Example3: Multi-block PCA with sparsity.
#' out <- moss(sim_blocks[-4],
#'      method="mbpca",
#'      dg.grid.right = seq(1,200,by=10),
#'      dg.grid.left = seq(1,100,by=2),
#'      alpha.right = 0.5,
#'      alpha.left = 1,
#'      tSNE=TRUE,
#'      clus=TRUE,
#'      clus.lab=lab.sub,
#'      plot=TRUE)
#'  out$clus.obj
#'  #This shows the 'weight' each omic block has on the variability explained by each PC.
#'  # Weights in each PC add up to one.
#'  out$block_weights
#'
#' #Example4: Partial least squares with sparsity (PLS).
#' out <- moss(sim_blocks[-4],
#'      K.X=500,
#'      K.Y=2,
#'      method="pls",
#'      dg.grid.right = seq(1,100,by=2),
#'      dg.grid.left = seq(1,100,by=2),
#'      alpha.right = 1,
#'      alpha.left = 1,
#'      tSNE=TRUE,
#'      clus=TRUE,
#'      clus.lab=lab.feat[1:2e3],
#'      resp.block=3,
#'      plot=TRUE,axes.pos=-(1:2))
#'  out$clus.obj
#'  table(out$sparse$u[,1] != 0,lab.feat[1:2000])
#'  table(out$sparse$v[,1] != 0,lab.feat[2001:3000])
#'
#' #Example5: PCA-LDA
#' out <- moss(sim_blocks,
#'      method="pca-lda",
#'      clus=TRUE,
#'      resp.block=4,
#'      clus.lab=lab.sub,
#'      plot=TRUE)
#' out$clus.obj
#'  }

moss <- function(data.blocks, scale.arg=TRUE, norm.arg=TRUE,method="pca",resp.block=NULL,
                   K.X=5,K.Y=K.X,verbose=TRUE,ncores=1,
                   dg.grid.left = NULL, dg.grid.right=NULL,
                   alpha.right=1,alpha.left=1,plot=FALSE,clus=FALSE,
                   clus.lab=NULL,tSNE=NULL,axes.pos=1:K.Y,approx.arg=FALSE,exact.dg=FALSE, use_fbm=FALSE) {
  
  #Inputs need to be a list of data matrices.
  if(!is.list(data.blocks)) stop("Input has to be a list with omic blocks.")
  
  #Only objects of class 'matrix' or 'FBM' are accepted.
  if(any(vapply(data.blocks, function (X) 
    any(vapply(c("matrix","FBM","array"), 
               function (x) inherits(X,x),TRUE)),TRUE) == FALSE))
  stop("Elements in data.blocks have to be 'array', 'matrix',or 'FBM' objects")
  
  #Number of data blocks.
  M <- length(data.blocks)
  
  if (is.null(tSNE) == FALSE | clus == TRUE | is.null(clus.lab) == FALSE) plot <- TRUE
  
  #Turning omic blocks into FBM.
  if (use_fbm == TRUE) {
    if (!requireNamespace("bigstatsr",quietly = TRUE)) stop("Package 'bigstatsr' needs to be installed to handle FBM objects.")
    else {
      if (verbose) message("Turn omic block into FBM objects.")
      data.blocks <- lapply(data.blocks, bigstatsr::as_FBM)
    }
  }
  
  #Checking the class of each data block.
  block.class <- rep("matrix", M)
  block.class[vapply(data.blocks, inherits, TRUE, what="FBM")] <- "FBM"
  
  #Printing method's name.
  c("Reduced Rank Regression (RRR)",
    "Partial Least Squares (PLS)",
    "Multi-Block Reduced Rank Regression (MBRRR)",
    "Multi-Block Partial Least Squares (MBPLS)",
    "Principal Components Analysis (PCA)",
    "Multi-Block Principal Components Analysis (MBPCA)",
    "Principal Components Analysis - Linear Discriminat Analysis (PCA-LDA)",
    "Multi-Block Principal Components Analysis - Linear Discriminant Analysis (MBPCA-LDA)",
    "Partial Least Squares - Linear Discriminant Analysis (PLS-LDA)",
    "Multi-Block Partial Least Squares - Linear Discriminant Analysis (MBPLS-LDA)",
    "Reduced Rank Regression - Linear Discriminant Analysis (RRR-LDA)",
    "Multi-Block Reduced Rank Regression - Linear Discriminant Analysis (MBRRR-LDA)"
  )[
    c('rrr','pls','mbrrr','mbpls',"pca","mbpca","pca-lda","mbpca-lda","pls-lda","mbpls-lda","rrr-lda","mbrrr-lda") == method] -> method.name
  
  title. <- paste0("| ",method.name," for ",M,ifelse(M == 1," data block and "," data blocks and "),K.Y,ifelse(M == 1," latent factor |"," latent factors |"))
  
  if (verbose) {
    message(rep("-", nchar(title.)))
    message(title.)
    message(rep("-", nchar(title.)))
  }
  
  #Checking if the right packages are present to handle approximated SVDs.
  if (approx.arg == TRUE) {
    if (any(block.class == "FBM"))  {if (!requireNamespace("bigstatsr",quietly = TRUE)) stop("Package 'bigstatsr' needs to be installed to handle FBM objects.")}
    else {if (!requireNamespace("irlba",quietly = TRUE)) stop("Package 'irlba' needs to be installed to get fast truncated SVD solutions.")}
  }

  #Checking if the right packages are present for Rtsne.
  if (is.null(tSNE) == FALSE) if(!requireNamespace("Rtsne",quietly = TRUE)) stop("Package 'Rtsne' needs to be installed to generate t-SNE projections.")

  #Checking if the right packages are present for dbscan.
  if (clus == TRUE) if(!requireNamespace("dbscan",quietly = TRUE)) stop("Package 'dbscan' needs to be installed for clustering.")

  #Checking if the right packages are present for plotting.
  if (plot == TRUE) {
    if(!requireNamespace("viridis",quietly = TRUE)) stop("Package 'viridis' needs to be installed for graphical displays.")
    if(!requireNamespace("ggplot2",quietly = TRUE)) stop("Package 'ggplot2' needs to be installed for graphical displays.")
    if(!requireNamespace("ggpmisc",quietly = TRUE)) stop("Package 'ggpmisc' needs to be installed for showing peaks on the PEV trajectory.")
  }
                             
  #Available methods.
  if (!any(method %in% c("pca","mbpca","pls","mbpls","rrr","mbrrr","pca-lda","mbpca-lda","pls-lda","mbpls-lda","rrr-lda","mbrrr-lda"))) {
    stop("Method ",method," not supported. Try one of these: pca, mbpca, pls, mbpls, rrr, mbrrr, pca-lda, mbpca-lda, pls-lda, mbpls-lda, rrr-lda, or mbrrr-lda.")
  }

  #At least two latent factor for tSNE.
  if (is.null(tSNE) == FALSE & K.Y < 2) stop("Number of latent factors needs to be larger or equal than 2 for tSNE projection.")

  #Only non-missing data accepted.
  if (verbose) message("Checking for missing values.")
  if (any(unlist(lapply(data.blocks, function(x) prepro_na(x) > 0)))) stop("Needs to imput NA's before runing R")

  #If method is a supervised one, the first response block is chosen: Y = data.blocks[[1]]
  if (!any(method %in% c("pca","mbpca"))) {
    if (is.null(resp.block)) {
      warning("If not specified, the first data block in the list will be used as response!")
      resp.block <- 1
    }
    else {
      #The response block needs to correspond to one element in the list
      if (!any(resp.block %in% 1:length(data.blocks))) stop("resp.block needs to be an integer between 1 and the total number of data blocks")
      else {data.blocks <- data.blocks[c(resp.block,(1:length(data.blocks))[-resp.block])];
      message("Block ",resp.block, " used as response.")}
    }
  }
  if (!grepl(method,pattern = "-lda")) if (any(lapply(data.blocks,function(x) inherits(x[,1],'numeric')) == FALSE)) stop("All blocks need to have numeric data.")

  #Naming data blocks if necessary.
  if (is.null(names(data.blocks))) names(data.blocks) <- paste("Block",1:M)
  else {
    tmp <- names(data.blocks) == ""
    names(data.blocks)[tmp] <- paste("Block",1:M)[tmp]
  }

  #Standardizing/Normalizing data blocks.
  if (scale.arg == T | norm.arg == T) {
    if (verbose) message("Standardizing/Normalizing data blocks.")
    if (grepl(method,pattern = "-lda")) data.blocks[-1] <- lapply(data.blocks[-1],prepro_sub,scale.arg=scale.arg,norm.arg=norm.arg)
    else data.blocks <-
      lapply(data.blocks,prepro_sub,scale.arg=scale.arg,norm.arg=norm.arg)
  }

  if (verbose) {
    #Correspondence between degree of sparsity and number of features only for LASSO penalties (alpha=1).
    if (alpha.right < 1) {
      if (alpha.right == 0) message('Ridge smoothing of right eigenvectors (no feature selection).')
      else message('Elastic net smoothing & selection of right eigenvectors. Number of features selected increases with degree of sparsity and alpha.right.')
    }
    if (alpha.left < 1){
      if (alpha.left == 0) message('Ridge smoothing of left eigenvectors (no', ifelse(grepl(method,pattern = "pca"),'subjects','responses') ,'selection).')
      else message('Elastic net smoothing & selection of left eigenvectors. Number of', ifelse(grepl(method,pattern = "pca"),'subjects','responses'), 'selected increases with degree of sparsity and alpha.left.')
    }
  }

  if (K.Y == 1 & plot == T) {warning("Plots will only be generated for more than ONE laten factor"); plot <- F}

  #Getting SVD from chosen method.
  if (method=="pca") {
    SVD_X <- rrr_pca(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD
    O <- SVD_X$Z
    left.lab <- "Subjects"
    right.lab <- "Predictors"
  }
  if (method=="mbpca") {
    SVD_X <- rrr_mb_pca(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$Z
    left.lab <- "Subjects"
    right.lab <- "Predictors"
  }
  if (method=="pls") {
    SVD_X <- rrr_pls(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="mbpls") {
    SVD_X <- rrr_mb_pls(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="rrr") {
    SVD_X <- rrr_rrr(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="mbrrr") {
    SVD_X <- rrr_mb_rrr(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="pca-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_lda(data.blocks,block.class,K.X,ncores,verbose,M)
    svd0 <- SVD_X$SVD
    O <- SVD_X$Z
    left.lab <- "Subjects"
    right.lab <- "Predictors"
  }
  if (method=="mbpca-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_mb_lda(data.blocks,block.class,K.X,ncores,verbose,M)
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$Z
    left.lab <- "Subjects"
    right.lab <- "Predictors"
  }
  if (method=="pls-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_pls(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    SVD_X$SVD$`w.x` <- SVD_X$SVD$u
    svd0 <- SVD_X$SVD
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="mbpls-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_mb_pls(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    SVD_X$SVD[[1]]$`w.x` <- SVD_X$SVD$u
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Response"
  }
  if (method=="rrr-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_rrr(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    SVD_X$SVD$`w.x` <- SVD_X$SVD$u
    svd0 <- SVD_X$SVD
    O <- SVD_X$LR
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }
  if (method=="mbrrr-lda") {
    if (ncol(data.blocks[[1]]) != 1 | inherits(data.blocks[[1]][,1], "character") == FALSE) stop("Response block needs to be a column with characters representing different classes.")
    if (min(table(data.blocks[[1]])) < 2) stop("Response block needs at least 2 samples by class")
    data.blocks[[1]] <- stats::model.matrix(~ -1 + data.blocks[[1]])
    SVD_X <- rrr_mb_rrr(data.blocks,block.class,K.X,K.Y,ncores,verbose,M)
    SVD_X$SVD[[1]]$`w.x` <- SVD_X$SVD$u
    svd0 <- SVD_X$SVD[[1]]
    O <- SVD_X$LR
    aux <- alpha.left
    alpha.left <- alpha.right
    alpha.right <- aux
    aux <- dg.grid.left
    dg.grid.left <- dg.grid.right
    dg.grid.right <- aux
    left.lab <- "Predictors"
    right.lab <- "Responses"
  }

  out <- NULL
  out$B <- O
  n <- nrow(O)
  out$dense <- SVD_X$SVD

  #Scree plot.
  if (plot & K.Y > 2) {
    aux.scree <- data.frame(y=c(svd0$d,c(0,diff(svd0$d),c(0,diff(c(0,diff(svd0$d)))))),
                            x=rep(1:K.Y,times=3),
                            type=rep(c("Scree plot","Eigen Values\n first derivative","Eigen Values\n second derivative"),each=K.Y))
    aux.scree$type <- factor(aux.scree$type,levels=c("Scree plot","Eigen Values\n first derivative","Eigen Values\n second derivative"),ordered = T)
    suppressWarnings(out$scree.plot <- ggplot2::ggplot(aux.scree, ggplot2::aes_string(x="x", y="y")) +
                       ggplot2::geom_line(col = "#21908C80") +
                       ggplot2::geom_point(col ="#44015480",pch=15)+
                       ggplot2::facet_wrap(.~type,scales = "free") +
                       ggplot2::scale_x_continuous("PC index")+
                       ggplot2::scale_y_continuous("Eigenvalues\n and derivatives")+
                       ggplot2::theme_minimal())
  }
  if (is.null(dg.grid.left) == F | is.null(dg.grid.right) == F) {
    #Sparsity constraints.
    if (verbose)       message("Imposing sparsity constraints.")
    aux.svd <- ssvdEN_sol_path(O = O,svd.0 = svd0,
                    scale= scale.arg,center = scale.arg,
                    dg.grid.right = dg.grid.right,dg.grid.left = dg.grid.left,
                    n.PC = K.Y,alpha.f = alpha.right,alpha.s = alpha.left,
                    plot = plot,approx = approx.arg,
                    verbose = verbose,left.lab = left.lab,right.lab = right.lab,exact.dg = exact.dg)

    out$sparse <- aux.svd$SVD
    if (plot) out$tun_dgSpar.plot <- aux.svd$plot

  }

  if (plot) {
    if (is.null(tSNE)) {
    if (is.null(clus.lab)) aux.name <- rep(left.lab,n)
    else aux.name <- clus.lab
    #Should we cluster left factors (e.g. first two PC's)?
    if (grepl(method,pattern = '-lda')) {
      if (clus) {
        if (verbose) message("Getting clusters via DBSCAN.")
        out$clus.obj <- tsne2clus(list(Y=scale(svd0$w.x[,1:2])),
                                  labels = aux.name,
                                  aest = aest.f(aux.name),
                                  eps_range = c(0,4),eps_res = 100,
                                  xlab = "LDF1",ylab="LDF2",clus=TRUE)
      }
      else {
        out$PC1_2.plot <- tsne2clus(list(Y=scale(svd0$w.x[,1:2])),
                                    labels = aux.name,
                                    aest = aest.f(aux.name),
                                    xlab = "LDF1",ylab="LDF2",clus=FALSE)
      }
    }
    else {
      if (clus) {
        if (verbose) message("Getting clusters via DBSCAN.")
        out$clus.obj <- tsne2clus(list(Y=scale(svd0$u[,1:2])),
                                  labels = aux.name,
                                  aest = aest.f(aux.name),
                                  eps_range = c(0,4),eps_res = 100,
                                  xlab = "PC1",ylab="PC2",clus=TRUE)
      }
      else {
       out$PC1_2.plot <- tsne2clus(list(Y=scale(svd0$u[,1:2])),
                                    labels = aux.name,
                                    aest = aest.f(aux.name),
                                    xlab = "PC1",ylab="PC2",clus=FALSE)
      }
    }
  }
    else {
    if (verbose) message("Calculating a tSNE map")
    if (grepl(method,pattern = '-lda')) {
      if (is.list(tSNE)) {tSNE$"Z" <- svd0$w.x; tSNE <- do.call(pca2tsne, tSNE)}
      else tSNE <- do.call(pca2tsne, list("Z"=svd0$w.x,"perp"=50,"n.samples"=1,"n.iter"=1e3))
      #Should we cluster left factors after tSNE?
      if (clus) {
        if (is.null(clus.lab)) aux.name <- rep(left.lab, n)
        else aux.name <- clus.lab
        if (verbose) message("Getting clusters via DBSCAN.")
        out$clus.obj <- tsne2clus(tSNE,
                                  labels = aux.name,
                                  aest = aest.f(aux.name),
                                  eps_range = c(0,4),eps_res = 100,
                                  xlab = paste0("tSNE_x{LDF",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),
                                  ylab=paste0("tSNE_y{LDF",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),clus=TRUE)
      }
      else {
        if (is.null(clus.lab)) aux.name <- rep(left.lab,n)
        else aux.name <- clus.lab
        out$tSNE.plot <- tsne2clus(tSNE,
                                   labels = aux.name,
                                   aest = aest.f(aux.name),
                                   xlab = paste0("tSNE_x{LDF",1,"-",ncol(data.blocks[[1]]),"}"),
                                   ylab=paste0("tSNE_y{LDF",1,"-",ncol(data.blocks[[1]]),"}"),clus = FALSE)


      }
    }
    else {
      if (is.list(tSNE)) {tSNE$"Z" <- svd0$u[, axes.pos]; tSNE <- do.call(pca2tsne, tSNE)}
      else tSNE <- do.call(pca2tsne, list("Z"=svd0$u[,axes.pos],"perp"=50,"n.samples"=1,"n.iter"=1e3))
      #Should we cluster left factors after tSNE?
      if (clus) {
        if (is.null(clus.lab)) aux.name <- rep(left.lab, n)
        else aux.name <- clus.lab
        if (verbose) message("Getting clusters via DBSCAN.")
        out$clus.obj <- tsne2clus(tSNE,
                                  labels = aux.name,
                                  aest = aest.f(aux.name),
                                  eps_range = c(0,4),eps_res = 100,
                                  xlab = paste0("tSNE_x{PC",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),
                                  ylab=paste0("tSNE_y{PC",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),clus=TRUE)
      }
      else {
        if (is.null(clus.lab)) aux.name <- rep(left.lab,n)
        else aux.name <- clus.lab
        out$tSNE.plot <- tsne2clus(tSNE,
                                   labels = aux.name,
                                   aest = aest.f(aux.name),
                                   xlab = paste0("tSNE_x{PC",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),
                                   ylab=paste0("tSNE_y{PC",paste0(range((1:K.Y)[axes.pos]),collapse ="-"),"}"),clus = FALSE)


      }
    }

    }
  }

  if(grepl(method,pattern = "mb")) {
    if (method == "mbpca") {M1 = M + 1;block.names1 <- names(data.blocks)}
    else {
      if (method == "mbpca-lda") {M1 = M; block.names1 <- names(data.blocks)[-1]}
      else {M1 = M;block.names1 <- names(data.blocks)[-1]}
    }
    names(out$dense) <- c("Global",block.names1)
    if(plot) {

      block.w <- data.frame(b=do.call("c", lapply(2:M1,
                                                  function(i) out$dense[[i]]$b[1:K.Y])),
                            Omic = rep(block.names1,each=K.Y),
                            LX.m = rep(1:K.Y,times=M1 - 1))

      out$block_weights <- ggplot2::ggplot(block.w, ggplot2::aes_string(x="LX.m", y="b",col="Omic",shape="Omic")) +
        ggplot2::scale_shape_manual(values=c(14:(13 + M1)))+
        ggplot2::scale_color_manual(values = viridis::viridis(M1 - 1))+
        ggplot2::geom_point(size=2.5)+
        ggplot2::scale_x_continuous("PC index")+
        ggplot2::scale_y_continuous("Weights")+
        ggplot2::geom_line()+
        ggplot2::theme_minimal()+
        ggplot2::theme(legend.position = "top")
    }

    pc.aux <- out$dense$Global$u
  }
  else  pc.aux <- out$dense$u

  #Evaluating the overlap between principal components and pre-defined groups of subjects.
  if (plot == TRUE & is.null(clus.lab) == FALSE) {
    message("Evaluating overlap between subjects selected and pre-established labels.")
    options(na.action="na.pass")
    Z <- stats::model.matrix(~-1 + clus.lab)
    colnames(Z) <- sort(unique(clus.lab[!is.na(clus.lab)]))

    #Creating a data.frame to store association between clusters and labels.
    clus.w <- do.call("rbind",lapply(1 : ncol(Z) , function(i) {
      if (is.null(dg.grid.left) == FALSE & alpha.left > 0) {
        x2.res <- apply(as.matrix(out$sparse$u),2,function(x) {
          suppressWarnings(test <- stats::chisq.test(table(x != 0,Z[,i])))
          return(c("Est"=test$statistic))
        })
      }

      else {
        x2.res <- apply(as.matrix(pc.aux),2,function(x) {
          test <- stats::kruskal.test(x,Z[,i])
          return(c("Est"=test$statistic))
        })
      }
      x2.res <- as.data.frame(x2.res)
      x2.res <- cbind(x2.res,1:K.Y,colnames(Z)[i])
      colnames(x2.res) <- c("Est","PC","Label")
      return(x2.res)
    }))

    if (is.null(dg.grid.left) == FALSE & alpha.left > 0) test.name <- "Chi-square statistics"
    else test.name <- "Kruskal-Wallis statistics"

    out$subLabels_vs_cluster <- ggplot2::ggplot(data=clus.w, ggplot2::aes_string(x="PC", y="Est",col="Label",shape="Label")) +
      ggplot2::scale_shape_manual(values=c(8:(8 + ncol(Z))))+
      ggplot2::scale_color_manual(values = viridis::viridis(ncol(Z)+1,option="A")[-(ncol(Z)+1)])+
      ggplot2::geom_point(size=2.5)+
      ggplot2::scale_x_continuous("PC index")+
      ggplot2::scale_y_continuous(test.name)+
      ggplot2::geom_line()+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.position = "top")

  }

  return(out)
}
