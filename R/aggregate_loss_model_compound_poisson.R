# Date: 27/07/16
#' Estimate a compound Poisson distribution.
#' @param poisson_mean Expected value of Poisson random variable.
#' @param lambda A vector of mle's of rate parameters from exponential mixture distribution.
#' @param p A vector of mle's of mixing proportions from exponential mixture distribution.
#' @param xmax A large integer representing a truncation point in the domain of exponential
#' mixture distribution.
#' @param N A large integer representing number of discretization points.
#' @param tol A number n such that P(Nmax > n) <= \code{tol} where Nmax has Poisson
#' with \code{possion_mean} distribution.
#' @return A list containing the support \code{x} of the convolution distribution and
#' an estimate \code{density} of the compound Poisson distribution.
compound_poisson <-
  function(poisson_mean,
           lambda,
           p,
           xmax,
           N = 2 ^ 12 ,
           tol = 0.005) {
    Nmax <- qpois(1 - tol, poisson_mean)
    f_S <- matrix(0, nrow = N, ncol = Nmax)
    for (k in 1:Nmax) {
      f_S[, k] <-
        dpois(k, poisson_mean) * m_fold_convolution(k, lambda = lambda, p = p, xmax)$convolution
    }
    eff_S <- apply(f_S, 1, sum)

    res <-
      list(x = m_fold_convolution(1, lambda = lambda, p = p, xmax)$x,
           density = eff_S)
    return(res)
  }
