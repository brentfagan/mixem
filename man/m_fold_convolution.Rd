% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_loss_model_m_fold_convolution.R
\name{m_fold_convolution}
\alias{m_fold_convolution}
\title{Estimate the m-fold convolution of a mixture of exponential distributions using FFT.}
\usage{
m_fold_convolution(m, lambda, p, xmax, N = 2^12)
}
\arguments{
\item{m}{An integer representing number of random variables in convolution.}

\item{lambda}{A vector of mle's of rate parameters from exponential mixture distribution.}

\item{p}{A vector of mle's of mixing proportions from exponential mixture distribution.}

\item{xmax}{A large integer representing a truncation point in domain of exponential mixture distribution.}

\item{N}{A large integer representing number of discretization points.}
}
\value{
A list containing the support \code{x} of the convolution distribution and
an estimate \code{convolution} of the m-fold convolution of the mixture distribution.
}
\description{
Estimate the m-fold convolution of a mixture of exponential distributions using FFT.
}
