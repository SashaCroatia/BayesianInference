#' Hamiltonian Monte Carlo for Logistic Regression with Normal or Data-Driven Priors
#' Fits a binary logistic regression model using Hamiltonian Monte Carlo (HMC) via the \pkg{nimble} and \pkg{nimbleHMC} packages, with normal priors based on previous coefficient samples.
#'
#' @param X A numeric matrix of predictors (n × M), where \code{n} is the number of observations and \code{M} is the number of predictors.
#' @param y A numeric or integer vector of binary responses (0/1), length \code{n}.
#' @param niter Integer, total number of MCMC iterations to run. Default is 5000.
#' @param nburnin Integer, number of burn-in iterations to discard. Default is 1000.
#' @param nchains Integer, number of chains to run. Default is 1.
#' @param inits Optional numeric vector of initial values for \code{beta}. If \code{NULL}, defaults to a zero vector.
#' @param KDE Logical, whether to use data-driven priors based on a set of coefficient samples.
#'   If \code{FALSE} (default), priors are \eqn{N(0, 100^2)}.
#'   If \code{TRUE}, priors are \eqn{N(\mathrm{mean}(C_j), \mathrm{sd}(C_j)^2)} for each coefficient \eqn{j}.
#' @param C A numeric matrix of coefficient samples (n × M) used only if \code{KDE = TRUE}.
#'
#' @return A list containing posterior samples
#'
#' @export
hmc_norm <- function(X, y, niter = 5000, nburnin = 1000, nchains = 1, inits = NULL, KDE = FALSE, C = NA) {
    n <- nrow(X)
    M <- ncol(X)

    # 1) Define model
    #---------------
    library(nimble)
    if (KDE == FALSE) {
        log_reg_code <- nimbleCode({
            # Piors
            for (j in 1:M) {
                beta[j] ~ dnorm(mean = 0, sd = 100)
            }

            # Likelihood
            for (i in 1:n) {
                logit(p[i]) <- inprod(X[i, 1:M], beta[1:M])
                y[i] ~ dbern(p[i])
            }
        })
    } else {
        log_reg_code <- nimbleCode({
            # Priors
            for (j in 1:M) {
                beta[j] ~ dnorm(mean = mean(C[1:n, j]), sd = sd(C[1:n, j]))
            }

            # Likelihood
            for (i in 1:n) {
                logit(p[i]) <- inprod(X[i, 1:M], beta[1:M])
                y[i] ~ dbern(p[i])
            }
        })
    }

    # 2) Prepare data, constants, init values
    #---------------
    hmc_consts <- list(n = n, M = M)
    hmc_data <- list(y = y, X = X, C = C)

    if (is.null(inits)) {
        hmc_inits <- list(beta = rep(0, M))
    } else {
        hmc_inits <- list(beta = inits)
    }

    # 3) Create model
    #---------------
    model <- nimbleModel(
        code = log_reg_code,
        constants = hmc_consts,
        data = hmc_data,
        inits = hmc_inits,
        calculate = FALSE,
        buildDerivs = TRUE
    )

    # 4) Compile model and mcmc to C++ for speed
    #---------------
    compiled_model <- compileNimble(model)

    # 5) Configure model
    #---------------
    library(nimbleHMC)
    model_config <- configureHMC(model)
    model_config$removeSamplers() # remove all
    model_config$addSampler(target = "beta")

    # 6) Build HMC algo
    #---------------
    HMC <- buildHMC(model = compiled_model, hmcConf = model_config)

    # 7) Compile HMC algo
    #---------------
    compiled_HMC <- compileNimble(HMC, project = model)

    # 8) Run HMC sampler
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


#' Prior-Recursive Hamiltonian Monte Carlo for Logistic Regression
#' Fits a binary logistic regression model using a prior-recursive approach with Hamiltonian Monte Carlo (HMC). The dataset is partitioned into \code{k} groups and regression coefficients from each group are used to construct data-driven priors for the next group.
#'
#' @param k Integer, number of partitions to split the data into.
#' @param X A numeric matrix of predictors (n × M), where \code{n} is the number of observations and \code{M} is the number of predictors.
#' @param y A numeric or integer vector of binary responses (0/1), length \code{n}.
#' @param niter Integer, total number of MCMC iterations to run for each partition. Default is 5000.
#' @param nburnin Integer, number of burn-in iterations to discard for each partition. Default is 1000.
#'
#' @return A matrix of posterior samples from the final partition, with columns corresponding to regression coefficients.
#'
#' @importFrom stats density approxfun
#' @export
pr_hmc_norm <- function(k, X, y, niter = 5000, nburnin = 1000) {
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
      
      C <- hmc_norm(parts_X[[j]], unlist(list(parts_y[[j]])), niter = niter, nburnin = nburnin, nchains = 1, inits = NULL, KDE = FALSE, C = NA)$samples
    
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

      C <- hmc_norm(parts_X[[j]], unlist(list(parts_y[[j]])), niter = niter, nburnin = nburnin, nchains = 1, inits = modes, KDE = TRUE, C = C)$samples

    }
  }
  return(C)
}