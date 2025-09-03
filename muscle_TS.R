# The code proposes the new MUSCLE4TS algorithm on R^d-valued heavy tailed times series

# In what follows, we use the relation between n, b and k : b = floor((n/k)^(1/phi))
# where n is the number of steps, b the block length, k the number of blocks that
# we consider extreme in the disjoint case, and phi a tunable parameter.
# prop is always the list of proportions of data that we will consider as extremes
# in the dijsoint case then we have approximately k = n*(prop^(phi/(phi-1)))

# The number of data in the consecutive case differ. We have the relation : 
# k_dis = k_cons * floor(n / bc) / (n-bc+1)

# It is recommanded to use phi = 2 (i.e b = sqrt(n/k) ) which corresponds to the 
# choice made in the paper by Buritica, Mikosch and Wintenberger

# We assume that the considered times series (X_t) has a tail index 1 and is stored 
# as a matrix d * n where d is such that X_t is included in R^d and n is the number 
# of steps


library(rlang)
source("muscle.r") # Assuming muscle.r contains necessary functions like clusters_vector and occ_vectors


## Preprocessing Function
#
# Assumes we have a R^d-valued time series Z_1,...,Z_n which is given as a d x n
# matrix Z.
#
# @param Z A d x n matrix representing the time series.
# @param b The size of the 'blocks' to be concatenated. Default is 1.
# @return A matrix of size (d x b) x ((n - b) + 1) where each column is the 
#         concatenation of consecutive vectors Z_j, ..., Z_(j + b - 1).
preprocess <- function(Z, b = 1) {
  # Handle the special case where d = 1, converting a vector to a 1-row matrix.
  if (is.matrix(Z) == FALSE) {
    Z <- t(as.matrix(Z))
  }
  
  # If block size is 1, no concatenation is needed.
  if (b == 1) {
    return(Z)
  }
  
  # Get dimensions of the input matrix.
  n <- ncol(Z)
  d <- nrow(Z)
  
  # Initialize the output matrix with the first 'b' columns of Z.
  # The size will be (d x b) x ((n - b) + 1).
  X <- Z[, 1:(n - b + 1)]
  
  # Loop to add subsequent vectors to the concatenated matrix.
  # For each i from 2 to b, it appends the (i-1)-th shifted version of Z.
  for (i in seq2(2, b)) {
    X <- rbind(X, Z[, i:(n - b + i)])
  }
  
  # Remove row names to avoid clutter and simplify subsequent processing.
  rownames(X) <- NULL
  
  return(X)
}

  
  ## Time Series Clustering Functions
  
  ### Main Wrapper for Time Series Clustering
  #
  # Applies the "muscle" algorithm for time series to a given matrix X. It
  # chooses the optimal parameters (k, b, number of directions) based on the
  # penalized log-likelihood.
  #
  # @param X A d x n matrix of time series data.
  # @param prop A sequence of proportions of data to consider as extremes.
# @param graph A boolean indicating whether to plot the results. Default is TRUE.
# @param phi A tunable parameter for block length calculation. Default is 2.
# @return A list containing the optimal cluster directions, the optimal proportion,
#         the optimal block length, the number of directions, and their weights.
muscle_clusters_TS <- function(X, prop = seq(0.005, 0.15, by = 0.005), graph = TRUE, phi = 2) {
  # Ensure X is a matrix, handling the d=1 case.
  if (is.matrix(X) == FALSE) {
    X <- t(as.matrix(X))
  }
  
  d <- nrow(X)
  n <- ncol(X)
  p <- length(prop)
  
  # Apply the main time series muscle algorithm.
  result <- muscle_TS(X, prop = prop, graph = graph, phi = phi)
  
  # Find the minimum penalized log-likelihood to select the optimal parameters.
  minimum <- which.min(result[[1]][, 5])
  
  # Extract the optimal number of directions, proportion, and matrix of directions.
  s_hat <- result[[1]][minimum, 4]
  prop_hat <- result[[1]][minimum, 1]
  M <- result[[2]][[minimum]]
  
  # Take only the 's_hat' most important directions.
  M <- as.matrix(M[, 1:s_hat])
  
  # Get the optimal block length 'b'.
  b <- result[[1]][minimum, 2]
  
  # Calculate the weights (proportions) of each direction.
  # The weights are stored in the last row of the M matrix.
  weights <- M[(d * b + 1), ] / sum(M[(d * b + 1), ])
  
  # Return the results as a list.
  result <- list(M, prop_hat, b, s_hat, weights)
  return(result)
}

  
  ### Time Series Muscle Algorithm
  #
  # Computes the KL (Kullback-Leibler) minimizer for different proportions of
  # extreme blocks and returns the output from `KL_minimizer_TS`. It also
  # optionally plots the results.
  #
  # @param X A d x n matrix of time series data.
  # @param prop A sequence of proportions of data to consider as extremes.
  # @param graph A boolean indicating whether to plot the results. Default is TRUE.
  # @param phi A tunable parameter for block length calculation. Default is 2.
# @return The output from `KL_minimizer_TS`.
muscle_TS <- function(X, prop = seq(0.005, 0.15, by = 0.005), graph = TRUE, phi = 2) {
  d <- nrow(X)
  n <- ncol(X)
  
  # Call the core function to find the KL minimizer for each proportion.
  result <- KL_minimizer_TS(X, prop, matrix_clusters = TRUE, phi = phi)
  
  # If `graph` is TRUE, generate plots.
  if (graph == TRUE) {
    # Plot the evolution of the penalized log-likelihood with proportion.
    plot(result[[1]][, 1], result[[1]][, 5], type = 'l', xlab = "proportion", 
         ylab = "penalized log-likelihood")
    
    # Plot the evolution of the number of extremal clusters with proportion.
    plot(result[[1]][, 1], result[[1]][, 4], type = 'l', xlab = "proportion", 
         ylab = expression(paste(hat(s)(k))))
  }
  
  return(result)
}

  
  ### Core KL Minimization Function
  #
  # Computes the penalized KL for different proportions of extreme blocks.
  #
  # @param X A d x n matrix of time series data.
  # @param prop A sequence of proportions of data to consider as extremes.
  # @param matrix_clusters A boolean indicating whether to return the list of
  #        direction matrices. Default is FALSE.
  # @param phi A tunable parameter for block length calculation. Default is 2.
  # @return A table with results for each proportion, and optionally a list of
#         matrices containing the directions and their occurrences.
KL_minimizer_TS <- function(X, prop = seq(0.005, 0.15, by = 0.005), matrix_clusters = FALSE, phi = 2) {
  n <- ncol(X)
  d <- nrow(X)
  length_prop <- length(prop)
  
  # Project consecutive blocks onto the positive l1 sphere for each proportion.
  mtx_clusters <- project_blocks_multiprop(X, prop, phi = phi)
  
  result <- c()
  
  # Loop through each proportion to compute the penalized KL.
  for (j in 1:length_prop) {
    # Ensure M is a matrix, handling the single direction case.
    M <- as.matrix(mtx_clusters[[j]])
    
    r <- ncol(M) # Number of unique directions.
    db <- nrow(M) - 1 # Value of d x b.
    occ_clusters <- M[nrow(M), ] # The number of occurrences for each direction.
    
    # Get current block length (bc), number of extremes (kc), and proportion (propc).
    bc <- (db / d)
    kc <- sum(occ_clusters)
    propc <- prop[j]
    nc <- ((n - bc) + 1) # Total number of consecutive blocks.
    
    # If there is only one direction, the KL is infinite.
    if (r == 1) {
      KL <- Inf
    } else {
      # Calculate the KL for different numbers of directions (s).
      # This is based on the log-likelihood of a multinomial distribution.
      KL <- 1:(r - 1) - lfactorial(kc) + kc * log(kc) + sum(lfactorial(occ_clusters)) -
        cumsum(occ_clusters[-r] * (log(occ_clusters)[-r])) - 
        (kc - cumsum(occ_clusters[-r])) * log((kc - cumsum(occ_clusters[-r])) / ((r - 1):1))
    }
    
    # Find the number of directions (`s_hat`) that minimizes the KL.
    s_hat <- which.min(KL)
    
    # Calculate the penalized KL, adding a penalty term for model complexity.
    KL_penalized <- (KL[s_hat]) / (kc) + kc / nc
    
    # Store the results.
    result <- rbind(result, c(propc, bc, r, s_hat, KL_penalized))
  }
  
  # Add column names for clarity.
  colnames(result) <- c("prop", "b", "r", "s", "KL penalized")
  
  # Return either the result table or a list including the cluster matrices.
  if (matrix_clusters == FALSE) {
    return(result)
  } else {
    return(list(result, mtx_clusters))
  }
}

  
  ### Projection of Blocks for Multiple Proportions
  #
  # For a sequence of proportions, this function projects consecutive blocks of a
  # time series onto the positive l1 sphere.
  #
  # @param X A d x n matrix of time series data.
  # @param prop A sequence of proportions.
  # @param phi A tunable parameter for block length calculation. Default is 2.
  # @return A list where each element is a matrix of projected unique directions
  #         and their occurrences for a given proportion.
project_blocks_multiprop <- function(X, prop, phi = 2) {
  # Ensure X is a matrix, handling the d=1 case.
  if (is.matrix(X) == FALSE) {
    X <- t(as.matrix(X))
  }
  
  n <- ncol(X)
  d <- nrow(X)
  
  # Calculate block lengths `b` for each proportion.
  b_seq <- (floor(prop**(-1 / phi)))
  
  result <- list()
  
  # Use `rle` (run-length encoding) to avoid redundant `preprocess` calls
  # when multiple proportions correspond to the same block length.
  rle_b <- rle(b_seq)
  length_cumsum = cumsum(rle_b$lengths)
  tmp <- 0
  
  for (i in 1:length(prop)) {
    # Get the current block length `bc`.
    bc <- rle_b$values[min(which(i <= length_cumsum))]
    
    # If the block length has changed, re-process the data into blocks.
    if (bc != tmp) {
      tmp <- bc
      Z <- preprocess(X, tmp)
    }
    
    # Determine the number of extremes `k` for the current proportion.
    k <- round(ncol(Z) * prop[i], 0)
    
    # Find the unique projected blocks for the current k.
    result[[i]] <- find_unique_projected_blocks(Z, k, d, bc)
  }
  
  return(result)
}

  
  ### Find Unique Projected Blocks
  #
  # Projects the top `k` most extreme blocks (based on l1 norm) onto the
  # positive l1 sphere and finds the unique directions.
  #
  # @param Z A matrix of consecutive blocks, as generated by `preprocess`.
  # @param k The number of most extreme blocks to consider.
  # @param d The dimensionality of the original data.
  # @param b The block length.
  # @return A matrix where each column is a unique projected direction. The last
#         row contains the number of occurrences for each direction.
find_unique_projected_blocks <- function(Z, k, d, b) {
  # Calculate the l1 norm of each block (column) of Z.
  Zl1 <- colSums(abs(Z))
  
  # Get the indices of the blocks sorted by their l1 norm in descending order.
  ord <- order(Zl1, decreasing = TRUE)
  
  # The radius of the sphere is the k-th largest l1 norm.
  u_k = Zl1[ord[k]]
  
  # Select the `k` most extreme blocks.
  Zps <- Z[, ord[(1:k)]]
  
  # Ensure Zps is a matrix, handling the case where k=1.
  if (is.matrix(Zps) == FALSE) {
    Zps <- t(Zps)
  }
  
  # Project each block onto the positive l1 sphere with radius `u_k`.
  Zps <- apply(Zps, 2, clusters_vector, c(u_k))
  
  # Ensure Zps is a matrix again after the `apply` call.
  if (is.matrix(Zps) == FALSE) {
    Zps <- t(Zps)
  }
  
  # Manually set columns with l1 norm equal to `u_k` to all ones, as per
  # the algorithm's definition of projection for these edge cases.
  tmp <- which(Zl1[ord[(1:k)]] == u_k, arr.ind = TRUE)
  Zps[, tmp] <- rep(1, d * b)
  
  # Remove any columns with NA values (which shouldn't happen but is a good practice).
  Zps <- Zps[, !is.na(Zps[1, ])]
  
  # Ensure Zps is a matrix after subsetting.
  if (is.matrix(Zps) == FALSE) {
    Zps <- t(Zps)
  }
  
  # Align the non-zero parts of the blocks to the beginning of the vector.
  # This is crucial for correctly identifying unique directions when working
  # with consecutive blocks and their inherent overlaps.
  if (nrow(Zps) > 1) {
    Zps <- apply(Zps, 2, align_non_zero_blocks, d, b)
  }
  
  # Ensure Zps is a matrix again.
  if (is.matrix(Zps) == FALSE) {
    Zps <- t(Zps)
  }
  
  # Find the unique vectors and count their occurrences.
  Zps <- occ_vectors(Zps)
  
  # Order the unique directions from most common to least common.
  ord <- order(Zps[d * b + 1, ], decreasing = TRUE)
  Zps <- Zps[, ord]
  
  return(Zps)
}

  
  ### Align Non-Zero Blocks
  #
  # Shifts a column vector of blocks so that its first non-zero block starts
  # at the beginning of the vector. This helps in identifying identical
  # directions regardless of their starting position within the concatenated block.
  #
  # @param x A single column vector representing a concatenated block.
  # @param d The dimensionality of the original data.
  # @param b The block length.
  # @return A new vector with the non-zero blocks shifted to the beginning.
align_non_zero_blocks <- function(x, d, b) {
  tmp <- as.matrix(rep(0, b * d))
  # Find indices of non-zero elements.
  ind <- which(x != 0, arr.ind = TRUE)
  
  # Calculate the number of d-dimensional zero vectors to remove from the start.
  j <- (min(ind) - 1) %/% d
  
  # Shift the non-zero part of the vector to the beginning.
  tmp[1:(d * (b - j))] <- x[(d * j + 1):(d * b)]
  
  return(tmp)
}