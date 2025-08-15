#' Hamiltonian Monte Carlo for Logistic Regression
#' Runs a Hamiltonian Monte Carlo (HMC) sampler for binary logistic regression using the \pkg{nimble} and \pkg{nimbleHMC} packages.
#'
#' @param X A numeric matrix of predictors (n × M), where n is the number of observations and M is the number of predictors.
#' @param y A numeric or integer vector of binary responses (0/1), length n.
#' @param niter Integer, total number of MCMC iterations to run.
#' @param nburnin Integer, number of burn-in iterations to discard.
#' @param nchains Integer, number of chains to run.
#' @param inits Optional named list of initial values. If \code{NULL}, defaults to zero for all \code{beta}.
#'
#' @return A list containing posterior samples and summary statistics as returned by \code{\link[nimble]{runMCMC}}.
#'
#' @export
hmc <- function(X, y, niter = 5000, nburnin = 1000, nchains = 1, inits = NULL) {
  #1) Define model
  #---------------
  library(nimble)
  log_reg_code <- nimbleCode({
    #Piors
    for (j in 1:M) {
      beta[j] ~ dnorm(mean = 0, sd = 100)
    }
    
    #Likelihood
    for (i in 1:n) {
      logit(p[i]) <- inprod(X[i, 1:M], beta[1:M])
      y[i] ~ dbern(p[i]) 
    }
  })
  
  #2) Prepare data, constants, init values
  #---------------
  hmc_consts <- list(n = nrow(X), M = ncol(X))
  
  hmc_data <- list(y = y, X = X)
  
  if(is.null(inits)){
    hmc_inits <- list(beta = rep(0, ncol(X)))
  } else {
    hmc_inits <- inits
  }
  
  #3) Create model
  #---------------
  model <- nimbleModel(
    code = log_reg_code,
    constants = hmc_consts,
    data = hmc_data,
    inits = hmc_inits,
    calculate = FALSE,
    buildDerivs = TRUE
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
