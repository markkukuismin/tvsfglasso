% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tvsfglasso.R
\name{tvsfglasso}
\alias{tvsfglasso}
\title{tvsfglasso}
\usage{
tvsfglasso(
  Y = NULL,
  N = NULL,
  use_cor = FALSE,
  pos = NULL,
  rep = FALSE,
  alpha = 0.1,
  K = 2,
  kernel = "gaussian",
  h = NULL
)
}
\arguments{
\item{Y}{list of time-series data.}

\item{N}{the number of time points.}

\item{use_cor}{compute the correlation matrix. Default \code{FALSE}.}

\item{pos}{time points where network estimates are computed.}

\item{rep}{does the data contain biological replicates. Default \code{FALSE}.}

\item{alpha}{scale-free glasso tuning parameter. Default \code{0.1}.}

\item{K}{the number of iterations (reweighting steps). Default \code{2}.}

\item{kernel}{the kernel used. Default \code{gaussian}.}

\item{h}{the bandwidth. Default \code{5.848/N^(1/3)}.}
}
\value{
List of sparse precision matrices, degrees of freedom, and smoothed covariance matrices.
}
\description{
Computes the time-varying scale-free glasso estimate.
}
\examples{

set.seed(1)

n <- 1
p <- 50
N <- 100

Data <- generate_tv_sf_data(p = p, n = n, N = N)

X <- Data$X

length(X)
dim(X[[1]])

Y <- matrix(0, p, N)

for(i in 1:N){
 
Y[, i] <- X[[i]]
 
}

dim(Y)

pos <- c(1, floor(N/2), N)

res <- tvsfglasso(Y = Y, N = N, pos = pos, rep = FALSE) 

n <- 5

Data <- generate_tv_sf_data(p = p, n = n, N = N)

Y <- Data$X

length(Y)
dim(Y[[1]])

res <- tvsfglasso(Y = Y, N = N, pos = pos, rep = TRUE) 

}
