stochOpt <- structure(function( ### Stochastic algorithm
### Stochastic algorithm for solving \eqn{P(H(theta, X))=c}, where \eqn{H} is a known function and \eqn{X} can be simulated
                               n, ### length of a trajectory
                               HX, ### a function of theta with values in \eqn{\{0,1\}} such that \eqn{HX(theta)=H(theta, X)}
                               c, ### target value for \eqn{P(H(theta, X))}
                               theta0=1/2,
                               p=3/4) {
  ii <- 1
  thetas <- NULL
  theta <- theta0
  thetas <- c(thetas, theta)
  while (ii < n) {
    quantity <- HX(theta)
    theta <- theta-(quantity-c)/ii^p
    thetas <- c(thetas, theta)
    ii <- ii+1
  }
  thetas
}, ex=function() {
  B <- 1e3
  alpha <- 0.2
  m <- 1e3
  tau <- function(alpha) qnorm(1-min(alpha, 1)*(1:m)/m) ## Simes
  ##  sim <- simulateGaussianNullsFromFactorModel(m, B=B, flavor="equi", rho=0.2)
  ## Y <- sim$Y
  HX <- function(lambda, alpha=0.2) {
    h <- 1
    P <- Matrix(1, m, length(h))
    rho <- 0.2
    Y <- simulateFactorModelNullsFromSingularValuesAndLoadings(m, h, P, rho=rho)$Y
    ts <- sort(Y, decreasing=TRUE)
    crit <- sum(ts>=tau(alpha*lambda))
    return(crit>0)
  }
  s0 <- stochOpt(B, HX, c=0.2, theta0=1)
  plot(s0)

  ## Compare with dichotomy
  flavor <- "equi-correlated"
  rho <- 0.2
  sim <- simulateGaussianNullsFromFactorModel(m, B=B, flavor=flavor, rho=rho)
  mat <- sim$Y
  str(mat)
  image(sim$Sigma)

  alpha <- 0.2
  maxSteps <- 1e3
  res <- getJointFWERThresholds(mat, refFamily="Simes", alpha, maxSteps=m, kMax=m)
  str(res)
})
