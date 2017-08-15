# Date: 27/07/16
# This file contains a function to estimate the m-fold convolution of a mixture of exponentials;
# assumes m iid random variables X and estimates the density of S = X_1 +...+ X_m. Also contains a
# function to estimate a compound Poisson distribution; assumes M iid random variables X and
# estimates the density of S = X_1 +...+ X_M where M has Poisson(M) distribution.
#' Estimate the m-fold convolution of a mixture of exponential distributions using FFT.
#' @param m An integer representing number of random variables in convolution.
#' @param lambda A vector of mle's of rate parameters from exponential mixture distribution.
#' @param p A vector of mle's of mixing proportions from exponential mixture distribution.
#' @param xmax A large integer representing a truncation point in domain of exponential mixture distribution.
#' @param N A large integer representing number of discretization points.
#' @return A list containing the support \code{x} of the convolution distribution and
#' an estimate \code{convolution} of the m-fold convolution of the mixture distribution.
m_fold_convolution <-
  function(m,
           lambda,
           p,
           xmax,
           N = 2 ^ 12) {
    dx <- xmax / N
    i <- 0:(N - 1)
    x <- i * dx
    du <- 2 * pi / (N * dx)
    b <- N / 2 * du
    a <- -b
    u <- a + i * du
    m_fourier_mix <-
      (sapply(1:length(u), function(j)
        sum(p * lambda / (lambda - 1i * u[j])))) ^ m
    convolution <-
      Re(du / (2 * pi) * exp(1i * a * x) * fft(m_fourier_mix))
    res <-
      list(x = x,
           convolution = convolution)
    return(res)

  }
