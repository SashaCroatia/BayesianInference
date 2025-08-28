#' Recursive Prior Proposal via Partitioned Metropolis Sampling
#' Runs a recursive prior proposal process by partitioning the dataset into \code{k} groups and sequentially fitting Metropolis-Hastings models. The first partition is sampled using an MVN prior, while subsequent partitions use random proposals from the posterior of the previous partition
#'
#' @param k Integer. Number of partitions to split the data into.
#' @param X Numeric matrix of predictors (rows = observations, columns = features).
#' @param y Numeric vector or matrix of responses (one column).
#' @param N Integer. Number of Metropolis samples for the first partition (default = 20000).
#' @param M Integer. Number of predictors (default = \code{ncol(X)}).
#' @param a_target Numeric. Target acceptance rate for adaptive tuning (default = 0.234).
#' @param b1 Numeric vector of initial proposal means (default = \code{rep(0, M)}).
#' @param batch_schedule Frequency of batching dataset for MH
#'
#' @return
#' A numeric matrix \code{C} containing posterior samples from the last partition's Metropolis-Hastings run.
#'
#' @importFrom stats density approxfun
#' @export
prior_proposal <- function(k, X, y, N = 20000, M = ncol(X), a_target = 0.234, b1 = rep(0, M), batch_schedule = seq(10, N, by = 10)) {
  #Partition the data
  #Number of rows
  n <- nrow(X)
  
  # Get row group assignments
  part_indices <- cut(seq_len(n), breaks = k, labels = FALSE)
  
  # Split row indices instead of splitting X directly
  index_groups <- split(seq_len(n), part_indices)
  
  # Turn each group of indices into a matrix
  parts_X <- lapply(index_groups, function(idx) X[idx, , drop = FALSE])
  y <- matrix(y, ncol = 1)
  parts_y <- lapply(index_groups, function(idx) y[idx, , drop = FALSE])
  
  #Prior recursive implementation
  for (j in seq(1,k)) { 
    if (j == 1) {
      C <- metropolis(parts_X[[j]], parts_y[[j]], N = N, prior = TRUE, update = FALSE, KDE = FALSE, C = NA, den_funs = NA, a_target = a_target, b1 = b1, batch = TRUE, batch_schedule = batch_schedule)
    
    } else {
      
      #Update number of samples
      N <- dim(C)[1] #number of samples
      
      #Create KDE functions and MAP (modes) vals
      den_funs <- list()
      modes <- list()
      d <- 1 #iterative
      while (d <= dim(C)[2]) {
        kde <- density(C[,d])
        den_funs[[d]] <- approxfun(kde$x, kde$y) #approx fun
        modes[[d]] <- kde$x[which.max(kde$y)]
        d <- d+1
      }
      modes <- unlist(modes)
      
      C <- metropolis(parts_X[[j]], parts_y[[j]], N = N, prior = FALSE, update = TRUE, KDE = FALSE, C = C, den_funs = NA, a_target = a_target, b1 = modes, batch = TRUE, batch_schedule = batch_schedule)
      
    }
  }
  return(C)
}