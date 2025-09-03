# The code proposes the muscle algorithm on R^d
# "MUltivariate Sparse CLustering for Extremes"
# by Nicolas Meyer and Olivier Wintenberger
# It was originally written by the same two 
# academics and slightly amended to be able to
# reuse some parts of it in Muscle4TS

# rm(list=ls())
library(rlang)

################################################################################
# This function takes a vector v with nonnegative coordinates as argument
# and gives a vector listing the ordered thresholds for which v has a different
# number of nonnull coordinates
################################################################################
cutoff_l1_norm <- function(v){
  d <- length(v)
  sort_v <- sort(v, decreasing = TRUE)
  cumsum_sort_v <- cumsum(sort_v)
  cutoff <- cumsum_sort_v - (1:d) * c(sort_v[seq2(2,d)],0) # 'seq2' to deal with the case d=1
  # the thresholds = v_(1) + ... + v_(k) - k v_(k+1)
  # for the last coordinate, it is v_1 + ... + v_d = cumsum_sort_v[d]
  #                                                = l1 norm of v
  # the command c(sort_v[2:d],0) just add a factor 0 at the end of this vector
  return(cutoff)
  # cutoff[i] contains the largest radius so that the positive l1 sphere
  # is such that the projection of v onto it has exactly i positive coordinates
}


################################################################################
# This function takes a vector v and a vector of thresholds 'thresholds' and 
# gives the clusters which contains the mass of pi(v).
# The library rlang is used here for 'seq2'.
# Idea: start with the largest threshold, look at which coordinates are positive
# after the projection, then look at the second largest threshold, and replace 
# some '0' by '1' if the associated coordinate is now positive. And so on
################################################################################
clusters_vector <- function(v, thresholds){
  length_v <- length(v)
  abs_v <- abs(v)
  order_abs_v <- order(abs_v)
  cutoff_v <- cutoff_l1_norm(abs_v)
  # we apply 'cutoff' to the |v|
  
  length_thresholds <- length(thresholds)
  thres_sort <- sort(thresholds, decreasing = TRUE)

  relevant_thres <- which(thresholds < cutoff_v[length_v], arr.ind = TRUE)
  # the coordinates of the relevant thresholds i.e smaller than the l1 norm of v
  # (the other ones are too large to be used for the projection)
  length_relevant_thres <- length(relevant_thres)
  # the one that are too big provide a vector of NA
  # we do not look at them
  
  mtx_clusters <- matrix(rep(NA, length_v*length_thresholds), nrow=length_v)
  mtx_clusters[, relevant_thres] <- 1
  # a matrix with NA for the thresholds that are to large and 1 for the others
  
  for (j in relevant_thres){
    # this loop allows one to replace some 1 by 0
    nb_pos_coord <- sum(cutoff_v < thresholds[j]) + 1
    nb_null_coord <- length_v - nb_pos_coord
    mtx_clusters[order_abs_v[seq2(1,nb_null_coord)], j] <- 0 # the coordinates for which pi(v)=0
    # the k coordinates set to 0 by the projection are the k smallest ones of v
    
    # we need to use seq2 since there might be some sequences of the form 1:0 
    # for which we do not want to do anything
  }

  # we add the sign of the coordinates and we reorder them following the order of 'thresholds'
  mtx_clusters <- sign(v)*mtx_clusters
  return(mtx_clusters)
}


################################################################################
# This function takes a matrix and gives the number of occurrence of each column
################################################################################
occ_vectors <- function(mtx){
  nb_col <- ncol(mtx)
  occ <- seq(nb_col)
  
  for (i in seq2(1,nb_col-1)) {
    if (!(occ[i] < i)){ 
      # if the current column have already been associated with a previous one 
      # then there is nothing to do
      dup <- colSums( abs( mtx[ ,seq(i+1,nb_col),drop=F] - mtx[ ,i] ) ) == 0
      # the current column is compared with all the following columns
      # the test consists in computing the difference of the vectors
      # and then checking if the absolute value of this difference is equal to 0
      occ[which(dup)+i] <- occ[i] 
      # occ contains for each column the column index of its first appearance
      }
  }
  # 2 cas utile seulement si nrow(mtx) ==1 possible
  if (nrow(mtx)<2){
    result <- rbind(t(mtx[ ,unique(occ)]), table(occ))
  } else {
    result <- rbind(as.matrix(mtx[ ,unique(occ)]), table(occ))
  }
  # only one copy of each different column is kept
  # the last line corresponds to the number of occurence
  # so the matrix 'result" is (d+1)*unique(occ)
  colnames(result) <- NULL
  return(result)
}


################################################################################
# This function takes a matrix X and a vector of thresholds 'thresholds'
# and gives for each threshold a matrix corresponding to
# 'clusters_vector' with the occurrence of each column
################################################################################
clusters_matrix <- function(X, thresholds){
  n <- ncol(X)
  d <- nrow(X)
  length_thresholds <- length(thresholds)
  clusters_all_thresholds <- apply(X, 2, clusters_vector, thresholds)
  # clusters_all_thresholds[(j-1)*d + 1 ):(j*d),i] contains the projection
  # of X[i] for the j-th threshold (i.e thresholds[j])
  if (is.matrix(clusters_all_thresholds)==FALSE){
    clusters_all_thresholds <- t(clusters_all_thresholds)
  } # to deal with the one dimensional case
  result <- list()
  
  for (j in 1:length_thresholds){
    mtx_clusters <- clusters_all_thresholds[( (j-1)*d + 1 ):(j*d) , !is.na(clusters_all_thresholds[(j-1)*d+1, ]) , drop=FALSE]
    # we remove the columns of NA of the considered block
    mtx_clusters <- occ_vectors(mtx_clusters)
    ordered_clusters <- order(mtx_clusters[d+1, ], decreasing=TRUE)
    mtx_clusters <- as.matrix(mtx_clusters[ ,ordered_clusters]) 
    # we sort the clusters in decreasing number of appearances
    result[[j]] <- mtx_clusters
  }
  result
}


################################################################################
# This function gives the evolution of the optimal value s_hat and k_hat
# for a data matrix X and a vector of proportion
################################################################################
KL_minimizer <- function(X, thresholds, matrix_clusters=FALSE){
  # This function gives for each k the value of s_hat which minimizes the KL
  # and the value of this minimizer
  # matrix_clusters = TRUE if we want the list of matrix M
  n <- ncol(X)
  d <- nrow(X)
  length_thresholds <- length(thresholds)
  mtx_clusters <- clusters_matrix(X, thresholds)
  
  result <- c()
  for (j in 1:length_thresholds){
    k <- sum(colSums(abs(X)) > thresholds[j])
    M <- mtx_clusters[[j]]
    r <- ncol(M)
    occ_clusters <- M[(d+1),]
    if (r==1){
      KL <- Inf
    }else{
      KL <- 1:(r-1) - lfactorial(k) + k*log(k) + sum(lfactorial(occ_clusters)) -
      cumsum(occ_clusters[-r]*(log(occ_clusters)[-r])) - (k-cumsum(occ_clusters[-r]))*log( (k-cumsum(occ_clusters[-r])) / ((r-1):1) )
    }
    # we compute the KL for different values of s
    s_hat <- which.min(KL)
    KL_penalized <- (KL[s_hat])/k + k/n
    result <- rbind(result, c(k, r, s_hat, KL_penalized))
  }
  colnames(result) <- c("k", "r", "hat_s", "KL penalized")
  if (matrix_clusters==FALSE){
    return(result)}else{
      return(list(result,mtx_clusters))
    }
}
  
  
################################################################################
# This function gives the evolution of the optimal value s_hat and k_hat
# for a data matrix X and a vector of proportion
################################################################################
muscle <- function(X, prop = seq(0.005,0.15,by=0.005), graph=TRUE){
  d <- nrow(X)
  n <- ncol(X)
  p <- length(prop)
  
  l1_norm <- apply(abs(X), 2, sum)
  sort_l1_norm <- sort(l1_norm, decreasing = TRUE)
  thresholds <- sort_l1_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*Proportion
  
  result <- KL_minimizer(X, thresholds)
  if (graph==TRUE){
    # Evolution of the KL with k
    plot(result[ ,1], result[ ,4], type='l', xlab = "k", ylab = "penalized log-likelihood")
    # Evolution of the number of extremal clusters with k
    plot(result[ ,1], result[ ,3], type='l', xlab = "k", ylab = expression(paste(hat(s)(k))))
  }
  return(result)
}


##################################################
# This function gives the extremal clusters for
# the optimal value of s_hat et hat_k
##################################################
muscle_clusters <- function(X, prop = seq(0.005,0.15,by=0.005)){
  d <- nrow(X)
  n <- ncol(X)
  p <- length(prop)
  
  l1_norm <- apply(abs(X), 2, sum)
  sort_l1_norm <- sort(l1_norm, decreasing = TRUE)
  thresholds <- sort_l1_norm[round(n*prop+1, 0)] # this corresponds to the list of thresholds u_n
  # we use the function 'round' since otherwise there are some approximations in n*prop
  
  result <- KL_minimizer(X, thresholds, matrix_clusters = TRUE)
  
  # Then we choose the minimal value for KL_penalised
  # and the associated values for M, s_hat, k_hat, etc.
  minimum <- which.min(result[[1]][,4])
  s_hat <- result[[1]][minimum,3]
  k_hat <- result[[1]][minimum,1]
  M <- result[[2]][[minimum]]
  M <- as.matrix(M[, 1:s_hat]) # we take only the s_hat first columns
  u <- sort_l1_norm[round(k_hat + 1, 0)] # the 'optimal' threshold
  weights <- M[(d+1), ]/ sum(M[(d+1), ])
  return(list(M, k_hat, u, s_hat, weights))
}

