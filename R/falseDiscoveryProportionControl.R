##' falseDiscoveryProportionControl
##'
##' False Discovery Proportion control
##'
##'
##' @param pval A vector of \eqn{p}-values
##' @param thr A non-increasing family of thresholds (on the \eqn{p}-value
##' scale) that controls joint FWER at level \eqn{(1-alpha)}
##' @param gamma A value in [0,1], the target FDP
##' @return A value in \eqn{t[hat] \in [0,1]} such that
##' \eqn{P(FDP(t[hat])>\gamma) \leq 1-alpha}.
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @examples
##'
##' m <- 200
##' alpha <- 0.1
##' thr <- alpha*1:m/m  ## Simes
##'
##' x <- rnorm(m)
##' pi0 <- .8
##' m1 <- round(m*(1-pi0))
##' cols <- 1+((1:m)<=m1)
##' idxs1 <- seq(from=1, to=m1, by=1)
##' mu <- 2
##' x[idxs1] <- x[idxs1] + mu
##' pval <- 1-pnorm(x)
##' gamma <- 0.2
##'
##' FDPbar <- falseDiscoveryProportionEnvelope(pval, thr)
##' tHat <- falseDiscoveryProportionControl(pval, thr, gamma)
##'
##' curve(FDPbar(x), 0, 0.005, type="S", ylim=c(0,1))
##' abline(h=gamma, col="lightgray")
##' stripchart(pval, add=TRUE, at=0, col=cols)
##' abline(v=thr, col="brown", lty=2, lwd=0.2)
##' abline(v=tHat, col="blue")
##'
falseDiscoveryProportionControl <- structure(function(
### False Discovery Proportion control
                                                       pval,
### A vector of \eqn{p}-values
                                                       thr,
### A non-increasing family of thresholds (on the \eqn{p}-value
### scale) that controls joint FWER at level \eqn{(1-alpha)}
                                                       gamma
### A value in [0,1], the target FDP
                                                       ) {
  ## kHat <- max(which(Vbar(thr) <= gamma*R(thr)))
  FDPbar <- falseDiscoveryProportionEnvelope(pval, thr)
  kHat <- max(which(FDPbar(thr)<=gamma))
  tHat <- thr[kHat]
### A value in \eqn{t[hat] \in [0,1]} such that \eqn{P(FDP(t[hat])>\gamma) \leq 1-alpha}.
}, ex=function(){
  m <- 200
  alpha <- 0.1
  thr <- alpha*1:m/m  ## Simes

  x <- rnorm(m)
  pi0 <- .8
  m1 <- round(m*(1-pi0))
  cols <- 1+((1:m)<=m1)
  idxs1 <- seq(from=1, to=m1, by=1)
  mu <- 2
  x[idxs1] <- x[idxs1] + mu
  pval <- 1-pnorm(x)
  gamma <- 0.2

  FDPbar <- falseDiscoveryProportionEnvelope(pval, thr)
  tHat <- falseDiscoveryProportionControl(pval, thr, gamma)

  curve(FDPbar(x), 0, 0.005, type="S", ylim=c(0,1))
  abline(h=gamma, col="lightgray")
  stripchart(pval, add=TRUE, at=0, col=cols)
  abline(v=thr, col="brown", lty=2, lwd=0.2)
  abline(v=tHat, col="blue")
})

############################################################################
## HISTORY:
## 2013-04-09
## o Created.
############################################################################

