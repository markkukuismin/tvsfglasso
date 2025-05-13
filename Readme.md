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

## Example 1

Here simulated data are generated with 100 variables, 200 time points, and 5 replicates,

```r
set.seed(1)

p <- 100
N <- 200
n <- 5

Data <- generate_tv_sf_data(p = p, n = n, N = N)

Y <- Data$X

pos <- 1:N

S_t <- kernel_cov_rep(X = Y, N = N, pos = pos)$S_t

Theta_est_1 <- sf_glasso(S = S_t[[1]], K = 2)

A1 <- Theta_est_1

A1[A1 != 0] <- 1

G1 <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected", diag = FALSE)

lo <- igraph::layout_with_fr(G1)

plot(G1, vertex.label = NA, layout = lo)

Theta_est_N <- sf_glasso(S = S_t[[N]], K = 2)

AN <- Theta_est_N

AN[AN != 0] <- 1

GN <- igraph::graph_from_adjacency_matrix(AN, mode = "undirected", diag = FALSE)

plot(GN, vertex.label = NA, layout = lo)
```