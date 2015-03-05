library(mvtnorm)
rho <- .2
sigma <- diag(2)
sigma[2, 1] <- rho
sigma[1, 2] <- rho
pmvnorm(mean=rep(0, 2), sigma, lower=rep(-Inf, 2), upper = c(1, 1))
