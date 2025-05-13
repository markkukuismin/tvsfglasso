#' @title sf_glasso
#'
#' @description Computes scale-free glasso estimate
#'
#' @param S list of (smoothed) covariance matrices.
#' @param alpha the adaptive tuning parameter. Default \code{0.1}.
#' @param K the number of iterations (reweighting steps). Default \code{2}.
#'
#' @return A scale-free glasso estimate.
#'
#' @examples
#'
#' library(glassoFast)
#' library(igraph)
#' library(MASS)
#'
#' n <- 5
#' p <- 50
#' N <- 100
#' 
#' Data <- generate_tv_sf_data(p = p, n = n, N = N)
#'
#' Y <- Data$X
#'
#' pos <- 1:N
#'  
#' S_t <- kernel_cov_rep(X = Y, N = N, pos = pos)$S_t
#' 
#' Theta_est <- sf_glasso(S = S_t[[1]], K = 2)
#'
#' @export
#' @importFrom glassoFast glassoFast
#' 
sf_glasso <- function(S = NULL, alpha = 0.1, K = 2){
  
  if(is.null(S)) stop("S is missing with no default")
  
  k <- 1
  
  p <- ncol(S)
  
  Theta_n <- diag(1, p)
  
  while(k <= K){
    
    ep <- diag(Theta_n)
    
    row_n = rowSums(abs(Theta_n) - diag(diag(abs(Theta_n))))
    
    xr <- 1/(row_n + ep)
    lambdaM <- outer(xr, xr, FUN = "+")
    lambdaM <- alpha*lambdaM
    
    beta <- 2*alpha/ep
    
    diag(lambdaM) <- beta
    
    results <- glassoFast::glassoFast(S, rho = lambdaM)
    
    Theta_n <- results$wi
    
    k <- k + 1
    
  }
  
  Theta_n
  
}