# Date: 08/01/17
# This file contains functions necessary to compute the maximum likelihood estimates of a mixture
# of exponential distributions when the data is grouped
#' Compute the maximum likelihood estimates of an exponential mixture distribution and number of mixtures when
#' the data is grouped.
#' @param x_grp Grouped data vector.
#' @param boundaries Interval boundaries of the grouped data.
#' @param data_type Three possible \code{data_type}s are \code{'interval', 'cross tab', or 'tabular'}
#' @param prec Precision used for testing convergence of mle's as measured by
#' squared difference of loglikelihood functions.
#' @param prec2 Estimated mixing proportions below \code{prec2} will cause algorithm to stop.
#' @param itermax Maximum number of EM algorithm iterations.
#' @param rmax Maximum number of mixtures to find.
#' @return A list of mixing proportions \code{p} and associated means \code{mean}.
em_group_fit <-
  function(x_grp,
           boundaries = NULL,
           data_type,
           prec = 1e-08,
           prec2 = 1e-04,
           itermax = 100,
           rmax = 10) {
    data <- format_data(x_grp,data_type)
    df <- data.frame(data[1])
    boundaries <- unlist(data[2])
    flag <- 1
    r <- 1 #initialize number of mixtures
    xmax <- boundaries[length(boundaries)]
    xmin <- (boundaries[1] + boundaries[2]) / 2
    for (i in 1:itermax) {
      if (r == 1) {
        initial_p <- rep(1 / r, r)
        initial_lambda <- seq(1 / xmax, 1 / xmin, length.out = r)
        iter <- 0
        diff <- 1
        lambda <- initial_lambda
        p <- initial_p
        for (j in 1:itermax) {
          old_lambda <- lambda
          old_p <- p
          lambda <-
            grouped_lambda_update(old_p, old_lambda, df, boundaries)
          p <- grouped_p_update(old_p, old_lambda, df, boundaries)
          diff <-
            (
              grouped_loglike(old_p, old_lambda, df, boundaries) - grouped_loglike(p, lambda, df, boundaries)
            ) ^ 2
          if (diff < prec) {
            break
          }
        }
        r <- r + 1
        new_lambda <- lambda
        new_p <- p
      }
      else{
        initial_p <- rep(1 / r, r)
        initial_lambda <- seq(1 / xmax, 1 / xmin, length.out = r)
        iter <- 0
        diff <- 1
        lambda <- initial_lambda
        p <- initial_p

        for (j in 1:itermax) {
          old_lambda <- lambda
          old_p <- p
          lambda <-
            grouped_lambda_update(old_p, old_lambda, df, boundaries)
          p <- grouped_p_update(old_p, old_lambda, df, boundaries)
          diff <-
            (
              grouped_loglike(old_p, old_lambda, df, boundaries) - grouped_loglike(p, lambda, df, boundaries)
            ) ^ 2
          if (diff < prec) {
            break
          }
        }
        candidate_p <- new_p
        candidate_lambda <- new_lambda
        flag <-
          grouped_loglike(p, lambda, df, boundaries) - grouped_loglike(new_p, new_lambda, df, boundaries)
        new_lambda <- lambda
        new_p <- p
        r <- r + 1
        if (flag < 0 || r > rmax || any(p < prec2)) {
          break
        }
      }
    }
    res <- list(p = candidate_p,
                lambda = candidate_lambda,
                means = 1 / candidate_lambda)
    return(res)
  }

## Mixture distribution function
grouped_mix_cdf <-
  function(p, lambda, x) {
    y <- p * (1-exp(-lambda %*% t(x)))
    l <- apply(y, 2, sum)
    return(l)
  }

## Grouped data log likelihood
grouped_loglike <- function(p,lambda,df,boundaries) {
  n_grp <- df[,ncol(df)]
  a <- boundaries[-length(boundaries)]
  b <- boundaries[-1]
  res <- n_grp*log(grouped_mix_cdf(p,lambda,b)-grouped_mix_cdf(p,lambda,a))
  return(sum(res))
}

## E_1 and E_2 are expectations required in EM algo
E_1 <- function(p,lambda,boundaries){
  a <- boundaries[-length(boundaries)]
  b <- boundaries[-1]
  res <- t(t(p*(exp(-lambda%*%t(a))-exp(-lambda%*%t(b))))/(grouped_mix_cdf(p,lambda,b)-grouped_mix_cdf(p,lambda,a)))
  return(res)
}

E_2 <- function(p,lambda,boundaries){
  a <- boundaries[-length(boundaries)]
  b <- boundaries[-1]
  res <- t(t(p/lambda*(exp(-lambda%*%t(a))*(lambda%*%t(a)+1)-
                         exp(-lambda%*%t(b))*(lambda%*%t(b)+1)))/(grouped_mix_cdf(p,lambda,b)-grouped_mix_cdf(p,lambda,a)))
  return(res)
}

## update functions
grouped_p_update <- function(p,lambda,df,boundaries){
  n <- sum(df[,ncol(df)])
  n_grp <- df[,ncol(df)]
  res <- rowSums(t(n_grp*t(E_1(p,lambda,boundaries))))/n
  return(res)
}

grouped_lambda_update <- function(p,lambda,df,boundaries){
  n <- sum(df[,ncol(df)])
  n_grp <- df[,ncol(df)]
  top <- rowSums(t(n_grp*t(E_1(p,lambda,boundaries))))
  bottom <- rowSums(t(n_grp*t(E_2(p,lambda,boundaries))))
  return(top/bottom)
}

## data format function
format_data <- function(x_grp,data_type){
  if (data_type == 'interval') {
    df <- data.frame(table(x_grp))
    boundary_list <-
      strsplit(as.character(gsub("\\(|\\]", "", df$x_grp)), ",")
    df$lower <- as.numeric(unlist(lapply(boundary_list, '[[', 1)))
    df$upper <- as.numeric(unlist(lapply(boundary_list, '[[', 2)))
    boundaries <- df$lower
    boundaries[nrow(df) + 1] <- df$upper[nrow(df)]
    return(list(df = df[,2],boundaries = boundaries))
  }
  else if (data_type == 'cross tab'){
    df <- data.frame(x_grp)
    boundary_list <-
      strsplit(as.character(gsub("\\(|\\]", "", df$x_grp)), ",")
    df$lower <- as.numeric(unlist(lapply(boundary_list, '[[', 1)))
    df$upper <- as.numeric(unlist(lapply(boundary_list, '[[', 2)))
    boundaries <- df$lower
    boundaries[nrow(df) + 1] <- df$upper[nrow(df)]
    return(list(df = df[,2],boundaries = boundaries))
  }
  else if (data_type == 'tabular'){ #need to build exception for case when upper limit = Inf
    df <- data.frame(x_grp)
    df$lower <- df[,1]
    df$upper <- df[,2]
    boundaries <- df$lower
    boundaries[nrow(df) + 1] <- df$upper[nrow(df)]
    return(list(df = df[,ncol(df)], boundaries = boundaries))
  }
  else{
    stop("Choose correct 'data_type'")
  }
}
