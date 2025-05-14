#' @title tvsfglasso
#'
#' @description Computes the time-varying scale-free glasso estimate.
#'
#' @param Y list of time-series data. 
#' @param N the number of time points.
#' @param use_cor compute the correlation matrix. Default \code{FALSE}.
#' @param pos time points where network estimates are computed.
#' @param rep does the data contain biological replicates. Default \code{FALSE}.
#' @param alpha scale-free glasso tuning parameter. Default \code{0.1}.
#' @param kernel the kernel used. Default \code{gaussian}.
#' @param K the number of iterations (reweighting steps). Default \code{2}.
#' @param h the bandwidth. Default \code{5.848/N^(1/3)}.
#'
#' @return List of sparse precision matrices, degrees of freedom, and smoothed covariance matrices.
#'
#' @examples
#' 
#' set.seed(1)
#'
#' n <- 1
#' p <- 50
#' N <- 100
#'
#' Data <- generate_tv_sf_data(p = p, n = n, N = N)
#'
#' X <- Data$X
#'
#' length(X)
#' dim(X[[1]])
#'
#' Y <- matrix(0, p, N)
#'
#' for(i in 1:N){
#'  
#' Y[, i] <- X[[i]]
#'  
#' }
#' 
#' dim(Y)
#'
#' pos <- c(1, floor(N/2), N)
#'
#' res <- tvsfglasso(Y = Y, N = N, pos = pos, rep = FALSE) 
#' 
#' n <- 5
#' 
#' Data <- generate_tv_sf_data(p = p, n = n, N = N)
#'
#' Y <- Data$X
#'
#' length(Y)
#' dim(Y[[1]])
#'
#' res <- tvsfglasso(Y = Y, N = N, pos = pos, rep = TRUE) 
#' 
#' @export
#' @importFrom stats cov dnorm sd cov2cor
#' @importFrom glassoFast glassoFast

tvsfglasso <- function(Y = NULL, N = NULL, use_cor = FALSE, pos = NULL, rep = FALSE, alpha = 0.1, K = 2, kernel = "gaussian", h = NULL){
  
  if(is.null(Y)) stop("Y is missing with no default")
  if(is.null(N)) stop("N is missing with no default")
  if(is.null(pos)) stop("pos is missing with no default")
  if(is.null(rep)) stop("Specify rep")
  
  if(is.null(h)) h <- 5.848/N^(1/3)
  
  if(!rep){
    
    S_t <- kernel_cov(X = Y, pos = pos, h = h, use_cor = use_cor, kernel = kernel)$S_t
    
    S_t <- lapply(seq(dim(S_t)[3]), function(x) S_t[ , , x])
    
  }
  
  if(rep){
    
    S_t <- kernel_cov_rep(X = Y, N = N, h = h, pos = pos, kernel = kernel)$S_t
    
  }
  
  p <- ncol(S_t[[1]])
  
  tvsfgl_res <- list(icov = list(), df = c(), St = list())
  
  for(i in 1:length(pos)){
    
    S_temp <- cov2cor(S_t[[i]])
    
    tvsfgl_res$icov[[i]] <- sf_glasso(S = S_temp, alpha = alpha, K = K)
    tvsfgl_res$df[i] = (sum(tvsfgl_res$icov[[i]] != 0) - ncol(tvsfgl_res$icov[[i]]))/2
    
  }
  
  tvsfgl_res$St <- S_t
  
  tvsfgl_res
  
}