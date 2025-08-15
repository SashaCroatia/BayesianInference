#' Custom Kernel Density Estimation Density Function for nimble
#' Provides a user-defined probability density function based on precomputed KDEs, for use within \pkg{nimble} models.
#'
#' @param x Numeric scalar, the value at which to evaluate the density.
#' @param iterative Numeric scalar (unused in current implementation, placeholder for model iteration indexing).
#' @param xkde Numeric vector of grid points (x-values) from the KDE.
#' @param ykde Numeric vector of density values (y-values) from the KDE corresponding to \code{xkde}.
#' @param log Logical (0 or 1), whether to return the log-density. Default is 0 (FALSE).
#'
#' @return A numeric scalar giving the density (or log-density if \code{log = 1}) at \code{x}.
#'
#' @import nimble
#' @export
d_kde <- nimble::nimbleFunction(
  run = function(x = double(0), iterative = double(0), xkde = double(1), ykde = double(1), log = integer(0, default = 0)) {

    kde_x <- xkde
    kde_y <- ykde

    #out of bounds
    if (x < min(kde_x) | x > max(kde_x)) return(1e-12)

    #find index
    i <- max(which(kde_x <= x))
    if (i == length(kde_x)) return(kde_y[i])

    #approxfun
    x0 <- kde_x[i]
    x1 <- kde_x[i+1]
    y0 <- kde_y[i]
    y1 <- kde_y[i+1]

    ans <- log(y0 + (y1 - y0)*(x - x0)/(x1 - x0))
    ans = x
    
    if(log) return(ans)
    return(exp(ans))
    returnType(double(0))
  },
  buildDerivs=TRUE
)


#' Custom KDE-Based Random Number Generator for nimble
#' Provides a user-defined random generation function intended for sampling
#' from a KDE within \pkg{nimble} models.
#'
#' @param n Integer scalar, the number of random samples to generate.
#' @param iterative Numeric scalar (unused in current implementation, placeholder for model iteration indexing).
#' @param xkde Numeric vector of grid points (x-values) from the KDE.
#' @param ykde Numeric vector of density values (y-values) from the KDE corresponding to \code{xkde}.
#'
#' @return A numeric scalar giving a random draw (currently uniform in \eqn{[0,1]}) as a placeholder.
#'
#' @import nimble
#' @export
r_kde<- nimble::nimbleFunction(
  run = function(n = integer(0), iterative = double(0), xkde = double(1), ykde = double(1)) {
    dev <- runif(1, 0, 1)
    return(dev)
    returnType(double(0))
  }
)

#' Hamiltonian Monte Carlo for Logistic Regression with Optional KDE Priors
#' Runs a Hamiltonian Monte Carlo (HMC) sampler for logistic regression  using \pkg{nimble} and \pkg{nimbleHMC}. Priors can be standard normal or custom kernel density estimates (KDE) based on posterior samples from a previous model fit.
#'
#' @param X Numeric matrix of predictors, dimension \eqn{n \times M}, where \eqn{n} is the number of observations and \eqn{M} is the number of predictors.
#' @param y Numeric vector of binary responses of length \eqn{n}.
#' @param niter Integer. Total number of MCMC iterations per chain (default = 5000).
#' @param nburnin Integer. Number of burn-in iterations to discard (default = 1000).
#' @param nchains Integer. Number of MCMC chains to run (default = 1).
#' @param inits Optional numeric vector of initial values for \code{beta}. If \code{NULL} (default), all coefficients are initialized at zero.
#' @param KDE Logical. If \code{FALSE} (default), standard normal priors are used for regression coefficients. If \code{TRUE}, KDE-based priors are constructed from \code{C}.
#' @param C Optional numeric matrix of posterior samples from a previous fit, used only when \code{KDE = TRUE}. Each column corresponds to a regression coefficient, and each row is a sample.
#'
#' @return A list containing:
#' \item{summary}{Summary statistics of posterior samples from all chains.}
#' \item{samples}{(If \code{summary = FALSE} in \code{\link{runMCMC}}) the raw samples.}
#'
#' @importFrom stats density
#' @export
hmc_kde <- function(X, y, niter = 5000, nburnin = 1000, nchains = 1, inits = NULL, KDE = FALSE, C = NA) {
  #0) Define constants
  #---------------
  n = nrow(X)
  M = ncol(X)
  
  #1) Define model
  #---------------
  library(nimble)
  if (KDE == FALSE) {
    log_reg_code <- nimbleCode({
      #Priors
      for (j in 1:M) {
        beta[j] ~ dnorm(mean = 0, sd = 100)
      }
      
      #Likelihood
      for (i in 1:n) {
        logit(p[i]) <- inprod(X[i, 1:M], beta[1:M])
        y[i] ~ dbern(p[i]) 
      }
    })
  } else {
    #Kde function lists
    n_grid <- 512
    kde_x <- matrix(NA, nrow = n_grid, ncol = M)
    kde_y <- matrix(NA, nrow = n_grid, ncol = M)

    for (j in 1:M) {
      kde <- density(C[, j], n = n_grid)
      kde_x[, j] <- kde$x  # each column is a KDE grid
      kde_y[, j] <- kde$y  # each column is the corresponding density
    }
    
    #Iterative list
    iter_list = c(1:M)
    x_kde_list = kde_x
    y_kde_list = kde_y
    

    log_reg_code <- nimbleCode({
      #Priors
      for (j in 1:M) {
        beta[j] ~ d_kde(iterative = iter[j], xkde = x_kde[1:n_grid, j], ykde = y_kde[1:n_grid, j])
      }
      
      #Likelihood
      for (i in 1:n) {
        logit(p[i]) <- inprod(X[i, 1:M], beta[1:M])
        y[i] ~ dbern(p[i]) 
      }
    })
  }
  
  #2) Prepare data, constants, init values
  #---------------
  if (KDE == FALSE) {
    hmc_consts <- list(n = n, M = M)
    hmc_data <- list(y = y, X = X)
    hmc_dim <- list(y = n, X = c(n, M))
  } else {
    hmc_consts <- list(n = n, M = M, iter = iter_list, n_grid = n_grid)
    hmc_data <- list(y = y, X = X, x_kde = x_kde_list, y_kde = y_kde_list)
    hmc_dim <- list(y = n, X = c(n,M), x_kde = c(n_grid,M), y_kde = c(n_grid,M))
  }
  
  if(is.null(inits)){
    hmc_inits <- list(beta = rep(0, M))
  } else {
    hmc_inits <- list(beta = inits)
  }
  
  #3) Create model
  #---------------
  model <- nimbleModel(
    code = log_reg_code,
    constants = hmc_consts,
    data = hmc_data,
    inits = hmc_inits,
    #calculate = FALSE,
    buildDerivs = TRUE,
    dimensions = hmc_dim
  )
  
  #4) Compile model and mcmc to C++ for speed
  #---------------
  compiled_model <- compileNimble(model)
  
  #5) Configure model
  #---------------
  library(nimbleHMC)
  model_config <- configureHMC(model)
  model_config$removeSamplers() #remove all
  model_config$addSampler(target = "beta")
  
  #6) Build HMC algo
  #---------------
  HMC <- buildHMC(model = compiled_model, hmcConf = model_config)
  
  #7) Compile HMC algo
  #---------------
  compiled_HMC <- compileNimble(HMC, project = model)
  
  #8) Run HMC sampler
  #---------------
  samples <- runMCMC(
    compiled_HMC,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    summary = TRUE
  )
  
  return(samples)
}


#' Partitioned Recursive HMC with KDE Priors
#' Runs a partitioned recursive Hamiltonian Monte Carlo (HMC) logistic regression, where the priors for later partitions are estimated using kdes from the posterior samples of earlier partitions.
#'
#' @param k Integer. Number of partitions to split the dataset into.
#' @param X Numeric matrix of predictors (rows = observations, columns = features).
#' @param y Numeric vector or single-column matrix of binary outcomes (0 or 1).
#' @param niter Integer. Total number of HMC iterations per partition (default = 5000).
#' @param nburnin Integer. Number of burn-in iterations to discard per partition (default = 1000).
#'
#' @return A numeric matrix of posterior samples from the last partition.  Columns correspond to regression coefficients (\code{beta}), and rows to iterations.
#'
#' @importFrom stats density
#' @export
pr_hmc_kde <- function(k, X, y, niter = 5000, nburnin = 1000) {
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
      
      C <- hmc_kde(parts_X[[j]], unlist(list(parts_y[[j]])), niter = niter, nburnin = nburnin, nchains = 1, inits = NULL, KDE = FALSE, C = NA)$samples
    
    } else {
      
      #Create KDE functions and MAP (modes) vals
      modes <- list()
      d <- 1 #iterative
      while (d <= dim(C)[2]) {
        kde <- density(C[,d])
        modes[[d]] <- kde$x[which.max(kde$y)]
        d <- d+1
      }
      modes <- unlist(modes)

      C <- hmc_kde(parts_X[[j]], unlist(list(parts_y[[j]])), niter = niter, nburnin = nburnin, nchains = 1, inits = modes, KDE = TRUE, C = C)$samples

    }
  }
  return(C)
}