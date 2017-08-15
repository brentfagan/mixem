# Date: 16/07/17
# This file contains functions necessary to compute the maximum likelihood estimates of a mixture
# of exponential distributions when the data is right censored
#' Compute the maximum likelihood estimates of an exponential mixture distribution and number of mixtures when
#' the data is right censored.
#' @param x A data vector.
#' @param cens_point Point at which data is censored.
#' @param n_c Number of right censored data points.
#' @param prec Precision used for testing convergence of mle's as measured by
#' squared difference of loglikelihood functions.
#' @param prec2 Estimated mixing proportions below \code{prec2} will cause algorithm to stop.
#' @param itermax Maximum number of EM algorithm iterations.
#' @param rmax Maximum number of mixtures to find.
#' @return A list of mixing proportions \code{p} and associated means \code{mean}.
em_censored_fit <- function(x,
                   cens_point,
                   n_c,
                   prec = 1e-08,
                   prec2 = 1e-04,
                   itermax = 100,
                   rmax = 10) {
  flag <- 1
  r <- 2 #initialize number of mixtures
  xmax <- cens_point
  xmin <- min(x)
  for(i in 1:itermax){
    if (r == 2) {
      initial_p <- rep(1 / r, r)
      initial_theta <- seq(xmin, xmax, length.out = r)
      iter <- 0
      diff <- 1
      theta <- initial_theta
      p <- initial_p

      for(j in 1: itermax){
        old_theta <- theta
        old_p <- p
        theta <- censored_theta_update(old_p, old_theta, x, cens_point, n_c)
        p <- censored_p_update(old_p, old_theta, x, cens_point, n_c)
        diff <- (loglike_cens(old_p, old_theta, x, n_c, cens_point) - loglike_cens(p, theta, x, n_c, cens_point)) ^ 2
        if(diff < prec){
          break
        }
      }
      r <- r + 1
      new_theta <- theta
      new_p <- p
    }
    else{
      initial_p <- rep(1 / r, r)
      initial_theta <- seq(xmin, xmax, length.out = r)
      iter <- 0
      diff <- 1
      theta <- initial_theta
      p <- initial_p

      for(j in 1: itermax){
        old_theta <- theta
        old_p <- p
        theta <- censored_theta_update(old_p, old_theta, x, cens_point, n_c)
        p <- censored_p_update(old_p, old_theta, x, cens_point, n_c)
        diff <- (loglike_cens(old_p, old_theta, x, n_c, cens_point) - loglike_cens(p, theta, x, n_c, cens_point)) ^ 2
        if(diff < prec){
          break
        }
      }
      candidate_p <- new_p
      candidate_theta <- new_theta
      flag <- loglike_cens(p, theta, x, n_c, cens_point) - loglike_cens(new_p, new_theta, x, n_c, cens_point)
      new_theta <- theta
      new_p <- p
      r <- r + 1
      if( flag < 0 || r > rmax || any(p < prec2)){
        break
      }
    }
  }
  res <- list(p = candidate_p, mean = candidate_theta)
  return(res)
}

##ungrouped data loglikelihood
loglike_scale <- function(p, theta, x) {
  v <- p * exp((-1/theta) %*% t(x))/theta
  inside <- apply(v, 2, sum)
  y <- sum(log(inside))
  return(y)
}

## Mixture distribution function
G <-
  function(p, theta, x) {
    y <- p * (1-exp((-1/theta) %*% t(x)))
    l <- apply(y, 2, sum)
    return(l)
  }

## censored data log likelihood
loglike_cens <- function(p,theta,x_k,n_cens,cens_point,right_cens=TRUE) {
  ungrp <- loglike_scale(p,theta,x_k)
  if (right_cens){
    res <- n_cens*log(G(p,theta,Inf)-G(p,theta,cens_point))
  }
  else{
    res <- n_cens*log(G(p,theta,cens_point)-G(p,theta,0))
  }
  return(sum(res) + ungrp)
}

## note that these have been calculated for right censored data
t_c <- function(p, theta, cens_point){
  y <- p*exp(-cens_point/theta)/sum(p*exp(-cens_point/theta))
  return(y)
}

t_x <- function(p, theta, x_k){
  y <- t(t(p/theta * exp((-1/theta) %*% t(x_k)))/colSums(p / theta * exp((-1/theta) %*% t(x_k))))
  return(y)
}

censored_p_update <- function(p, theta, x_k, cens_point, n_c){
  y <- (apply(t_x(p,theta,x_k), 1, sum) + n_c*t_c(p,theta,cens_point))/(length(x_k)+n_c)
  return(y)
}

censored_theta_update <- function(p, theta, x_k, cens_point, n_c){
  y <- (apply(t(t(t_x(p,theta,x_k))*x_k), 1, sum) + n_c*t_c(p,theta,cens_point)*(cens_point+theta))/
    (apply(t_x(p,theta,x_k), 1, sum) + n_c*t_c(p,theta,cens_point))
  return(y)
}
