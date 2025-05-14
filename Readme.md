# tvsfglasso: Time-varying scale-free graphical lasso

## Description

The R package `tvsfglasso` can be used to estimate dynamic networks from time-series data. See the example below.

For more details on the MDSD, please refer to: "tvsfglasso: time-varying scale-free graphical lasso (to appear). Supplementary data and code needed to reproduce the results reported in the article are available at: [https://github.com/markkukuismin/tvsfglasso_supplementary](https://github.com/markkukuismin/tvsfglasso_supplementary).

## Dependencies

Please make sure to install the following packages before using R package `tvsfglasso`. R with version later than 4.3.1 is needed.
```r
install.packages(c("igraph", "glassoFast"))
```

## Installation

The R package `tvsfglasso` can be installed from source files in the GitHub repository (R package `devtools` is needed):
```r
library(devtools)
install_github(repo="markkukuismin/tvsfglasso")
```

You can run tvglasso and tvsfglasso either by,

1. First computing the weighted covariance matrix and using them as input of glasso or scale-free glasso algorithms.

2. Use functions `tvglasso` and `tvsfglasso` of the `tvsfglasso` package.

The first option facilitates backtracking in the event of errors.

## Example 1: run tvsfglasso manually

Here simulated data are generated with 100 variables, 200 time points, and 5 replicates.

```r
set.seed(1)

p <- 100
N <- 200
n <- 5

Data <- generate_tv_sf_data(p = p, n = n, N = N)

Y <- Data$X
```

Estimate the time-varying network at each time-points. First, compute the weighted covariance matrix

```r
pos <- 1:N

S_t <- kernel_cov_rep(X = Y, N = N, pos = pos)$S_t
```

Next, compute the scale-free graphical lasso (glasso) estimate. The resulting precision matrices are interdependent and exhibit a scale-free structure.

Here estimates are computed at time-points 1 and N,

```r
S_temp <- cov2cor(S_t[[1]]) 

Theta_est_1 <- sf_glasso(S = S_temp, K = 2)
```

Alternatively you can use any other glasso estimator if you do not want the resulting precision matrix to exhibit a scale-free structure.

```r
#S_temp <- cov2cor(S_t[[1]])
#Theta_est_1 <- glassoFast::glassoFast(S = S_temp, rho = 0.1)
```

```r
A1 <- Theta_est_1

A1[A1 != 0] <- 1

G1 <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected", diag = FALSE)

lo <- igraph::layout_with_fr(G1)

plot(G1, vertex.label = NA, layout = lo)
```

```r
S_temp <- cov2cor(S_t[[N]]) 

Theta_est_N <- sf_glasso(S = S_temp, K = 2)

AN <- Theta_est_N

AN[AN != 0] <- 1

GN <- igraph::graph_from_adjacency_matrix(AN, mode = "undirected", diag = FALSE)

plot(GN, vertex.label = NA, layout = lo)
```

## Example 2: use tvsfglasso

Estimate the time-varying network at each time-point.

```r
pos <- 1:N

res <- tvsfglasso(Y = Y, N = N, pos = pos, rep = TRUE)
```

Alternatively you can use `tvglasso` if you do not want the resulting precision matrix to exhibit a scale-free structure.

```r
#res <- tvglasso(Y = Y, N = N, pos = pos, rep = TRUE)
```

The resulting adjacency matrices are identical to those computed earlier.

```r
A1_1 <- res$icov[[1]]

A1_1[A1_1 != 0] <- 1

all.equal(A1_1, A1)

G1_1 <- igraph::graph_from_adjacency_matrix(A1_1, mode = "undirected", diag = FALSE)

lo <- igraph::layout_with_fr(G1_1)

plot(G1_1, vertex.label = NA, layout = lo)
```

```r
AN_N <- res$icov[[N]]

AN_N[AN_N != 0] <- 1

all.equal(AN_N, AN)

GN_N <- igraph::graph_from_adjacency_matrix(AN_N, mode = "undirected", diag = FALSE)

plot(GN_N, vertex.label = NA, layout = lo)
```