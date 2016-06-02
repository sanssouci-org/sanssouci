##' falseDiscoveryProportionEnvelope
##'
##' confidence envelope for the False Discovery Proportion process
##'
##'
##' @param pval A vector of \eqn{p}-values
##' @param thr A non-increasing family of thresholds (on the \eqn{p}-value
##' scale) that controls joint FWER at level \eqn{(1-alpha)}
##' @return A function \eqn{f:[0,1] \mapsto [0,1]} giving a
##' \eqn{(1-alpha)}-confidence envelope for the FDP process.
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @export
##' @examples
##'
##' m <- 1000
##' alpha <- 0.1
##' thr0 <- alpha*1:m/m  ## Simes
##'
##' x <- rnorm(m)
##' pi0 <- .9
##' m1 <- round(m*(1-pi0))
##' cols <- 1+((1:m)<=m1)
##' idxs1 <- seq(from=1, to=m1, by=1)
##' mu <- 2
##' x[idxs1] <- x[idxs1] + mu
##' pval <- 1-pnorm(x)
##' gamma <- 0.2
##'
##' FDPbar <- falseDiscoveryProportionEnvelope(pval, thr0)
##' tHat <- falseDiscoveryProportionControl(pval, thr0, gamma)
##' t1 <- "FDP envelope using Simes' thresholds"
##' t2 <- sprintf("m=%s; alpha=%s", m, alpha)
##'
##' curve(FDPbar(x), 0, 5*tHat, type="S", ylim=c(0,1), ylab=expression(bar(FDP)(t)), xlab="t")
##' mtext(t1, side=3, adj=0)
##' mtext(t2, side=3, adj=1)
##' abline(h=gamma, col="lightgray")
##' stripchart(pval[idxs1], add=TRUE, at=0, col=2, cex=0.3)
##' stripchart(pval[-idxs1], add=TRUE, at=0, col=1, cex=0.3)
##' abline(v=thr0, col="brown", lty=2, lwd=0.2)
##' abline(v=tHat, col="blue")
##'
falseDiscoveryProportionEnvelope <- structure(function(
### confidence envelope for the False Discovery Proportion process
                                                       pval,
### A vector of \eqn{p}-values
                                                       thr
### A non-increasing family of thresholds (on the \eqn{p}-value
### scale) that controls joint FWER at level \eqn{(1-alpha)}
                                                       ) {
  ## sanity checks
  m <- length(pval)
  stopifnot(length(thr)==m)

  Vbar <- function(t) falseDiscoveryEnvelope(t, thr)
  R <- function(t) sapply(t, FUN=function(u) sum(pval<=u))
  FUN <- function(t) Vbar(t)/pmax(R(t), 1)
### A function \eqn{f:[0,1] \mapsto [0,1]} giving a \eqn{(1-alpha)}-confidence envelope for the FDP process.
}, ex=function(){
  m <- 1000
  alpha <- 0.1
  thr0 <- alpha*1:m/m  ## Simes

  x <- rnorm(m)
  pi0 <- .9
  m1 <- round(m*(1-pi0))
  cols <- 1+((1:m)<=m1)
  idxs1 <- seq(from=1, to=m1, by=1)
  mu <- 2
  x[idxs1] <- x[idxs1] + mu
  pval <- 1-pnorm(x)
  gamma <- 0.2

  FDPbar <- falseDiscoveryProportionEnvelope(pval, thr0)
  tHat <- falseDiscoveryProportionControl(pval, thr0, gamma)
  t1 <- "FDP envelope using Simes' thresholds"
  t2 <- sprintf("m=%s; alpha=%s", m, alpha)

  curve(FDPbar(x), 0, 5*tHat, type="S", ylim=c(0,1), ylab=expression(bar(FDP)(t)), xlab="t")
  mtext(t1, side=3, adj=0)
  mtext(t2, side=3, adj=1)
  abline(h=gamma, col="lightgray")
  stripchart(pval[idxs1], add=TRUE, at=0, col=2, cex=0.3)
  stripchart(pval[-idxs1], add=TRUE, at=0, col=1, cex=0.3)
  abline(v=thr0, col="brown", lty=2, lwd=0.2)
  abline(v=tHat, col="blue")
})

############################################################################
## HISTORY:
## 2013-04-09
## o Created.
############################################################################

