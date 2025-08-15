#' Recursive prior estimation with adaptive Metropolis sampling
#' Performs a recursive Bayesian update where data are split into \code{k} partitions. The posterior from one partition is used to form a kernel density estimate (KDE) prior for the next partition. Sampling for each stage is performed using the \code{\link{metropolis}} function.
#'
#' @param k Integer; number of data partitions.
#' @param X Numeric design matrix of predictors.
#' @param y Numeric response vector (binary outcomes 0/1).
#' @param N Integer; number of MCMC iterations per partition (default \code{20000}).
#' @param a_target Numeric; target acceptance rate for Metropolis sampler (default \code{0.234}).
#' @param b1 Numeric vector; initial parameter values for the first partition (default \code{rep(0, M)}).
#'
#' @details
#' The data are split row-wise into \code{k} roughly equal groups. For the first group, the Metropolis sampler is run with a MVN prior. For subsequent groups, KDEs are fit to the posterior samples from the previous stage, producing univariate priors for each parameter. The prior mean vector is set to the mode of each KDE for the starting values.
#'
#' @return Numeric matrix of posterior samples from the final stage after all recursive updates.
#'
#' @importFrom stats density approxfun
#' @export
prior_recursive <- function(k, X, y, N = 20000, M = ncol(X), a_target = 0.234, b1 = rep(0, M)) {
  #Partition the data
  #Number of rows and columns
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
      C <- metropolis(parts_X[[j]], parts_y[[j]], N = N, prior = TRUE, update = FALSE, KDE = FALSE, C = NA, den_funs = NA, a_target = a_target, b1 = b1, batch = TRUE)
    
    } else {
      
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
      
      C <- metropolis(parts_X[[j]], parts_y[[j]], N = N, prior = TRUE, update = FALSE, KDE = TRUE, C = C, den_funs = den_funs, a_target = a_target, b1 = modes, batch = TRUE)
      
    }
  }
  return(C)
}