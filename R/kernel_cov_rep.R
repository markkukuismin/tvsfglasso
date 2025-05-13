#' @title kernel_cov_rep
#'
#' @description Computes the smoothed covariance matrix estimates for biological replicates.
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
#' S_t <- kernel_cov_rep(X = Y, N = N, pos = pos)
#'
#' @export
#' @importFrom stats cov dnorm

kernel_cov_rep <- function(X = NULL, N = NULL, pos = NULL, h = NULL, kernel = "gaussian") {
  
  if(!is.list(X)) stop("X must be a list")
  if(is.null(X)) stop("X is missing with no default")
  if(is.null(N)) stop("N is missing with no default")
  
  p <- nrow(X[[1]])
  
  if(is.null(h)) h <- 5.848/N^(1/3)
  
  S_temp <- vector("list", N)
  
  S_t <- vector("list", length(pos))
  
  for(i in 1:N){
    
    S_temp[[i]] <- stats::cov(t(X[[i]]))
    
  }
  
  s = pos/max(pos)
  
  if(kernel == "epanechnikov") K <- function(u) (3/4)*(1 - u^2)
  if(kernel == "gaussian") K <- function(u) stats::dnorm(u)
  
  i <- 1
  
  for(tt in s){
    W <- K(abs(s - tt)/h)
    W <- W/sum(W)
    Sw <- matrix(0, p, p)
    j <- 1
    for(w in W){
      Sw <- Sw + w*S_temp[[j]]
      j <- j + 1
    }
    S_t[[i]] <- Sw
    i <- i + 1
  }
  
  result <- list(S_t = S_t)
  return(result)
}