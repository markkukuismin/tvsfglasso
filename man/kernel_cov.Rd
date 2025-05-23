% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_cov.R
\name{kernel_cov}
\alias{kernel_cov}
\title{kernel_cov}
\usage{
kernel_cov(
  X = NULL,
  pos = NULL,
  h = NULL,
  use_cor = FALSE,
  kernel = "gaussian"
)
}
\arguments{
\item{X}{the time-series data.}

\item{pos}{time points where network estimates are computed.}

\item{h}{the bandwidth. Default \code{5.848/N^(1/3)}.}

\item{use_cor}{compute the correlation matrix. Default \code{FALSE}.}

\item{kernel}{the kernel used. Default \code{gaussian}.}
}
\value{
List of smoothed covariance or correlation matrices.
}
\description{
Computes the smoothed covariance matrix estimates
}
\examples{

set.seed(1)

n <- 1
p <- 50
N <- 100

Data <- generate_tv_sf_data(p = p, n = n, N = N)

Y <- Data$X

length(Y)
dim(Y[[1]])

X <- matrix(0, p, N)

for(i in 1:N){
 
X[, i] <- Y[[i]]
 
}

pos <- 1:N

S_t <- kernel_cov(X = X, pos = pos)

}
