##' falseDiscoveryEnvelope
##'
##' confidence envelope for the False Discovery process
##'
##' If the input threshold family yields joint FWER less than \eqn{alpha}, then
##' the output (seen as a function of t) is a \eqn{(1-alpha)}-confidence
##' envelope for the False Discovery process.
##'
##' @param t A vector of values in [0,1]
##' @param thr A non-increasing family of thresholds (on the \eqn{p}-value
##' scale) that control joint FWER at level \eqn{(1-alpha)}
##' @return A vector of integers in [0, length(thr)].
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @examples
##'
##' m <- 200
##' alpha <- 0.2
##' thr <- alpha*1:m/m  ## Simes
##'
##' falseDiscoveryEnvelope(0.1, thr)
##' curve(falseDiscoveryEnvelope(x, thr))
##'
##' ## adapting to a specific dependency setting
##' B <- 1e3
##' rho <- 0.3
##' flavor <- c("independent", "equi-correlated", "3-factor")[2]
##' mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
##'
##' res <- jointFWERControl(mat, refFamily="Simes", alpha, maxSteps=100)
##' thr <- 1-pnorm(res$thr) ## converting to p-value scale
##' curve(falseDiscoveryEnvelope(x, thr), add=TRUE, col=2)
##'
falseDiscoveryEnvelope <- structure(function(
### confidence envelope for the False Discovery process
                                             t,
### A vector of values in [0,1]
                                             thr
### A non-increasing family of thresholds (on the \eqn{p}-value
### scale) that control joint FWER at level \eqn{(1-alpha)}
                                             ) {
  ##details<<If the input threshold family yields joint FWER less than
  ##\eqn{alpha}, then the output (seen as a function of t) is a
  ##\eqn{(1-alpha)}-confidence envelope for the False Discovery
  ##process.
  thr <- c(thr, 1)  ## implcitly capping FD to 1
  sapply(t, FUN=function(u) min(which(thr >= u)) -1)
### A vector of integers in [0, length(thr)].
}, ex=function(){
  m <- 200
  alpha <- 0.2
  thr <- alpha*1:m/m  ## Simes

  falseDiscoveryEnvelope(0.1, thr)
  curve(falseDiscoveryEnvelope(x, thr))

  ## adapting to a specific dependency setting
  B <- 1e3
  rho <- 0.3
  flavor <- c("independent", "equi-correlated", "3-factor")[2]
  mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)

  res <- jointFWERControl(mat, refFamily="Simes", alpha, maxSteps=100)
  thr <- 1-pnorm(res$thr) ## converting to p-value scale
  curve(falseDiscoveryEnvelope(x, thr), add=TRUE, col=2)
})

############################################################################
## HISTORY:
## 2013-04-05
## o Created.
############################################################################

