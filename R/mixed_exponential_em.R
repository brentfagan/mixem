# Date: 27/07/16
# This file contains functions necessary to compute the maximum likelihood estimates of a mixture
# of exponential distributions
#' Compute the maximum likelihood estimates of an exponential mixture distribution and number of mixtures.
#' @param x A data vector.
#' @param prec Precision used for testing convergence of mle's as measured by
#' squared difference of loglikelihood functions.
#' @param prec2 Estimated mixing proportions below \code{prec2} will cause algorithm to stop.
#' @param itermax Maximum number of EM algorithm iterations.
#' @param rmax Maximum number of mixtures to find.
#' @return A list of mixing proportions \code{p} and associated means \code{mean}.
em_fit <- function(x, prec = 1e-08, prec2 = 1e-04, itermax = 100, rmax = 10) {
  flag <- 1
  r <- 1 #initialize number of mixtures
  xmax <- max(x)
  xmin <- min(x)
  for(i in 1:itermax){
    if (r == 1) {
      initial_p <- rep(1 / r, r)
      initial_lambda <- seq(1/xmax, 1/xmin, length.out = r)
      iter <- 0
      diff <- 1
      lambda <- initial_lambda
      p <- initial_p

      for(j in 1: itermax){
        old_lambda <- lambda
        old_p <- p
        new_g <- g(old_p, old_lambda, x)
        new_f <- f(old_p, old_lambda, x)
        lambda <- lambda_update(new_f, new_g, x)
        p <- p_update(old_p, new_f, new_g)
        diff <- (loglike(old_p, old_lambda, x) - loglike(p, lambda, x)) ^ 2
        if(diff < prec){
          break
        }
      }
      r <- r + 1
      new_lambda <- lambda
      new_p <- p
    }
    else{
      initial_p <- rep(1 / r, r)
      initial_lambda <- seq(1/xmax, 1/xmin, length.out = r)
      iter <- 0
      diff <- 1
      lambda <- initial_lambda
      p <- initial_p

      for(j in 1: itermax){
        old_lambda <- lambda
        old_p <- p
        new_g <- g(old_p, old_lambda, x)
        new_f <- f(old_p, old_lambda, x)
        lambda <- lambda_update(new_f, new_g, x)
        p <- p_update(old_p, new_f, new_g)
        diff <- (loglike(old_p, old_lambda, x) - loglike(p, lambda, x)) ^ 2
        if(diff < prec){
          break
        }
      }
      candidate_p <- new_p
      candidate_lambda <- new_lambda
      flag <- loglike(p, lambda, x) - loglike(new_p, new_lambda, x)
      new_lambda <- lambda
      new_p <- p
      r <- r + 1
      if( flag < 0 || r > rmax || any(p < prec2)){
        break
      }
    }
  }
  res <- list(p = candidate_p, mean = candidate_lambda^(-1))
  return(res)
}
# Helper functions for main function em_fit
###########################################
# This function computes the loglikelihood of a mixture of exponentials
# Inputs:
#   p: a vector of mixing proportions
#   lambda: a vector of rate parameters for the exponential distribution
#   x: a data vector
# Outputs:
#   A scalar representing the value of the loglikelihood function
loglike <- function(p, lambda, x) {
  v <- p * lambda * exp(-lambda %*% t(x))
  inside <- apply(v, 2, sum)
  y <- sum(log(inside))
  return(y)
}

# This function computes the mixture density g
# Inputs:
#   p: a vector of mixing proportions
#   lambda: a vector of rate parameters for the exponential distribution
#   x: a data vector
# Outputs:
#   A vector g(x[1]), g(x[2]),..., g(x[n])
g <-
  function(p, lambda, x) {
    y <- p * lambda * exp(-lambda %*% t(x))
    l <- apply(y, 2, sum)
    return(l)
  }

# This function computes the exponential density f
# Inputs:
#   p: a vector of mixing proportions
#   lambda: a vector of rate parameters for the exponential distribution
#   x: a data vector
# Outputs:
#   A length(lambda) by length(x) matrix f_j(x[1]),...,f_j(x[n]) jth row
f <-
  function(p, lambda, x) {
    y <- lambda * exp(-lambda %*% t(x))
    return(y)
  }

# This function computes the updated mixing proportion in EM algorithm
# Inputs:
#   p: a vector of mixing proportions
#   lambda: a vector of rate parameters for the exponential distribution
#   x: a data vector
# Outputs:
#   A vector of updated mixing proportions
p_update <- function(p, f, g) {
  y <- t(f) / g
  l <- p * apply(y, 2, mean)
  return(l)
}

# This function computes the updated rate parameter in EM algorithm
# Inputs:
#   p: a vector of mixing proportions
#   lambda: a vector of rate parameters for the exponential distribution
#   x: a data vector
# Outputs:
#   A vector of updated rate parameters
lambda_update <- function(f, g, x) {
  y <- t(f) / g
  numerator <- apply(y, 2, sum)
  z <- x * t(f) / g
  denominator <- apply(z, 2, sum)
  l <- numerator / denominator
  return(l)
}
