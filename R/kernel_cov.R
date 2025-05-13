#' @title kernel_cov
#'
#' @description Computes the smoothed covariance matrix estimates
#'
#' @param X the time-series data. 
#' @param pos time points where network estimates are computed.
#' @param h the bandwidth. Default \code{5.848/N^(1/3)}.
#' @param use_cor compute the correlation matrix. Default \code{FALSE}.
#' @param kernel the kernel used. Default \code{gaussian}.
#'
#' @return List of smoothed covariance or correlation matrices.
#'
#' @examples
#'
#' library(huge)
#' library(igraph)
#' library(MASS)
#'
#' n <- 5
#' p <- 100
#' N <- 100
#' 
#' Data <- generate_tv_sf_data(p = p, n = n, N = N)
#'
#' Y <- Data$X
#'
#' pos <- 1:N
#'
#' S_t <- kernel_cov(X = Y, N = N, pos = pos)
#'
#' @export
#' @importFrom stats cov dnorm sd

kernel_cov <- function(X = NULL, pos = NULL, h = NULL, use_cor = FALSE, kernel = "gaussian") {
  
  if(!is.numeric(X)) stop("X must be numeric")
  if(is.null(X)) stop("X is missing with no default")
  
  p <- nrow(X)
  N <- ncol(X)
  
  if(is.null(pos)) pos <- 1:N
  
  if(is.null(h)) h <- 5.848/N^(1/3) # (200)^(1/3) = 5.848, See Zhou et al (2010) examples
  
  sdX <- rep(1, p)
  
  if(use_cor){
    for(i in 1:p){
      sdX[i] <- stats::sd(X[i, ])
      X[i, ] <- X[i, ]/sd(X[i, ])
    }
  }
  
  S_t <- array(0, c(p, p, N))
  
  for(i in 1:N){
    if(kernel == "epanechnikov") Kh <- pmax(3/4*(1 - ((pos - i)/((N - 1)*h))^2), 0)
    if(kernel == "gaussian") Kh <- pmax(dnorm(((pos - i)/((N - 1)*h))), 0) 
    omega <- Kh/sum(Kh)
    ind <- which(omega != 0)
    X_pos <- X[, pos[ind]]
    S_t[, , i] <- (X_pos*rep(omega[ind], rep(p, length(ind))))%*%t(X_pos)
  }
  
  result <- list(S_t = S_t, sdX = sdX)
  return(result)
}