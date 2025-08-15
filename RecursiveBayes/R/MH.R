#' Proposal distribution sampler
#' Draw a single sample from a multivariate normal distribution.
#'
#' @param mu Numeric vector; the mean vector.
#' @param lambda Numeric scalar; scaling factor for the covariance matrix.
#' @param sigma Numeric matrix; base covariance matrix.
#'
#' @return Numeric vector of sampled values.
#'
#' @export
proposal <- function(mu, lambda, sigma) {
  MASS::mvrnorm(1, mu = mu, Sigma = lambda * sigma)
}


#' Log-likelihood of the logistic regression function
#' Compute the log-likelihood of logistic regression parameters.
#'
#' @param beta Numeric vector of regression coefficients.
#' @param X Numeric design matrix.
#' @param y Numeric binary response vector (0/1).
#'
#' @return Numeric scalar: log-likelihood value.
#'
#' @export
logistic_ll <- function(beta, X, y) {
  z <- X %*% beta
  p <- 1 / (1 + exp(-z))
  loglike <- sum(log(p[y == 1])) + sum(log((1 - p)[y == 0]))
}


#' Log-prior density under a multivariate normal distribution
#' Compute the log-density for a given vector specified from an MVN proposal
#'
#' @param distro Numeric vector of parameter values.
#' @param mu Numeric vector; mean vector of the prior.
#' @param sigma Numeric covariance matrix of the prior.
#'
#' @return Numeric scalar: log-prior value
#'
#' @export
log_prior <- function(distro, mu, sigma) {
  prob <- mvtnorm::dmvnorm(distro, mean = mu, sigma = sigma)
  log(prob)
}


#' Log prior density using independent KDEs
#' Computes the log prior density for a parameter vector, assuming each parameter in the vector has an independent prior distribution estimated via kernel density estimation (KDE), one density function per parameter.
#'
#' @param beta Numeric vector of parameter values to evaluate.
#' @param C Numeric matrix of posterior samples or simulated values used to define the KDE support for each parameter. Each column corresponds to a parameter.
#' @param den_funs List of univariate density functions (e.g., from \code{\link[stats]{approxfun}} or \code{\link[stats]{density}}) giving the estimated prior density for each parameter. The list must be the same length as \code{beta}.
#'
#' @details
#' The function assumes parameter independence, so the joint prior density is the product of the univariate priors. If a parameter value lies outside the range of \code{C[, d]} for that parameter, the log prior is set to \code{-10000} (a large negative penalty) and no further terms are added.
#'
#' @return Numeric scalar: the sum of the log prior densities for all parameters, or a large negative value if any parameter is out of bounds.
#'
#' @export
log_prior_kde <- function(beta, C, den_funs) {
  prior_d = 0
  d <- 1 #iterative
  while (d <= dim(C)[2]) {
    if ( min(C[,d]) < beta[d] && beta[d] < max(C[,d]) ){
      density_at_point <- den_funs[[d]](beta[d])
      prior_d <- prior_d + log(density_at_point) #Assume independence
    }
    else {
      prior_d <- -10000
    }
    d <- d+1
  }
  return(prior_d)
}


#' Adaptive Metropolis-Hastings Sampler for Logistic Regression
#' Runs a (potentially adaptive) Metropolis-Hastings MCMC algorithm for sampling regression coefficients in a logistic regression model. Supports Gaussian priors or KDE-based empirical priors, and can update the proposal distribution adaptively in batches.
#'
#' @param X Numeric matrix of predictors (\eqn{n \times M}).
#' @param y Numeric vector or single-column matrix of binary responses (0/1).
#' @param N Integer. Total number of MCMC iterations (default = 20000).
#' @param prior Logical. If \code{TRUE}, includes a prior term in the posterior calculation (default = TRUE).
#' @param update Logical. If \code{TRUE}, proposed samples are taken from the provided chain \code{C} rather than the proposal distribution (default = FALSE).
#' @param KDE Logical. If \code{TRUE}, uses a KDE-based prior via \code{\link{log_prior_kde}}; otherwise uses Gaussian prior via \code{\link{log_prior}} (default = FALSE).
#' @param C Numeric matrix. Previous posterior samples used for KDE prior (if \code{KDE = TRUE}) or for replaying proposals if \code{update = TRUE}.
#' @param den_funs List of functions. Each function returns a density value for a given coefficient (only used if \code{KDE = TRUE}).
#' @param M Integer. Number of predictors (default = \code{ncol(X)}).
#' @param a_target Numeric. Target acceptance rate for adaptation (default = 0.234).
#' @param b1 Numeric vector of length \code{M}. Initial value of regression coefficients (default = \code{rep(0, M)}).
#' @param batch Logical. If \code{TRUE}, adapts proposal covariance in batches; otherwise uses single-step Robbins–Monro adaptation (default = TRUE).
#' @param batch_schedule Integer vector. Iterations at which batch updates occur (default = \code{seq(50, N, by = 50)}).
#'
#' @details
#' The algorithm follows an adaptive Metropolis-Hastings framework:
#' \enumerate{
#'   \item Propose a new sample from a multivariate normal proposal distribution (\code{\link{proposal}}).
#'   \item Compute the acceptance ratio using the logistic regression likelihood (\code{\link{logistic_ll}}) 
#'         plus the specified prior.
#'   \item Accept or reject the proposal according to the Metropolis criterion.
#'   \item Adapt the proposal distribution to improve mixing, either in batches or incrementally.
#' }
#' If \code{KDE = TRUE}, the prior is estimated from previous MCMC output \code{C} by fitting kernel density 
#' estimates for each coefficient separately.
#'
#' @return A numeric matrix of posterior samples for regression coefficients, with the first 10\% discarded as burn-in.
#'
#' @export
metropolis <- function(X, y, N = 20000, prior = TRUE, update = FALSE, KDE = FALSE, C = NA, den_funs = NA, M = ncol(X), a_target = 0.234, b1 = rep(0, M), batch = TRUE, batch_schedule = seq(50, N, by = 50)) {
  #0) Initial values
  if (batch == FALSE) {
    a_target <- a_target/2.26 
  }
  n_accept <- 0 #number of acceptances
  b <- matrix(0, nrow = N, ncol = M)
  b[1,] <- b1 #rep(0, M)
  
  #Keep history of proposals
  P_list <- rep(NA, N - 1)
  P_list[1] <- logistic_ll(b[1,], X, y)
  
  #Adaptation parameters
  batch_index <- 1
  Tb <- 0 #start of batch
  mu <- rep(0, M)
  Sigma <- diag(M)
  lambda <- (2.38**2)/M
  gam <- 1
  
  for (i in 2:N) {
    #1) Propose new position beta[i-1] -> beta[i]' from proposal distro
    if (update == FALSE) {
      b_prop <- proposal(mu, lambda, Sigma)
    } else {
      r <- sample(2:N, 1)
      b_prop <- C[r, ]
    }
    
    #2) Compute transition probability
    if (prior == TRUE) {
      if (KDE == FALSE) {
        P_list[i] <- logistic_ll(b_prop, X, y) + log_prior(b_prop, rep(0, M), diag(M))  #P(beta*)
      } else {
        P_list[i] <- logistic_ll(b_prop, X, y) + log_prior_kde(b_prop, C, den_funs)
      }
    } else {
      P_list[i] <- logistic_ll(b_prop, X, y)
    }
    
    #3) Acceptance ratio
    R <- exp(P_list[i] - P_list[i-1])
    alpha <- min(R, 1)
    
    #4) Generate a random number u_i from [0,1]. Accept or reject proposal
    u <- runif(1)  # u_i
    if (u <= alpha) {
        b[i,] <- b_prop  # Accept y as next state
        n_accept <- n_accept + 1
    } else {
        b[i,] <- b[i-1,]  # Stay at previous state of Markov chain
        P_list[i] <- P_list[i-1]
    }
    
    #5) Adaptive step
    if (batch == TRUE) {
      if (batch_index <= length(batch_schedule) && i == batch_schedule[batch_index]) {
        #Batch update
        Tb1 <- i #T_{b+1}
        b_batch <- b[(Tb + 1):Tb1, , drop = FALSE]
        n_batch <- nrow(b_batch)
        
        #Covar update
        cov_sum <- matrix(0, M, M)
        for (j in 1:n_batch) {
          diff <- matrix(b_batch[j,] - mu, ncol=1)
          cov_sum <- cov_sum + (diff%*%t(diff) - Sigma)
        }
        Sigma <- Sigma + cov_sum/n_batch + 1e-5*diag(M)
        
        #Mean update
        mean_diff <- colMeans(b_batch) - mu
        mu <- mu + mean_diff
        
        #Parameter update
        Tb <- Tb1
        batch_index <- batch_index + 1
        lambda <- exp( log(lambda) + (alpha-a_target)/n_batch )
      }
    } else { #Batch == FALSE
      gam <- 1/i
      Sigma <- Sigma + gam*( (b[i,]-mu)%*%t(b[i,]-mu) - Sigma )
      mu <- mu + gam*( b[i,]-mu )
      lambda <- exp( log(lambda) + gam*(alpha-a_target) )
    }
  }

  n <- as.integer(0.9 * N)
  C <- b[(N - n + 1):N,] #discard first 10% of samples
  
  cat("Acceptance rate:", n_accept/(N-1), "\n")
  return(C)
}