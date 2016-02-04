stepDownControl <- structure(function(
### get a family of thresholds that provide step down control of the joint FWER
    stat,
### A vector of \eqn{m} test statistics
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test statistics under the null hypothesis. \describe{
### \item{m}{is the number of tested hypotheses}
### \item{B}{is the number of Monte-Carlo samples}}
    refFamily=c("Simes", "kFWER", "LR06"),
### A character value or a function that returns a vector of \code{m} thresholds (see Details).
    H0=NULL,
### A numeric vector, the indices of true null hypotheses (between \code{1} and \code{m}) Defaults to \code{NULL}. This paramenter is typically used to perform simulations.
    ...,
### Further arguments to be passed to getJointFWERControl.
    verbose=FALSE
### If \code{TRUE}, print results of intermediate calculations.
### Defaults to \code{FALSE}.
    ) {
  m <- nrow(mat)
  stopifnot(length(stat)==m)
  stopifnot(refFamily=="kFWER")  ## other flavor not implemented yet (?)

  if (!is.null(H0)) {
      ## sanity checks
      stopifnot(length(H0)<=m)
      stopifnot(all(H0<=m))
  }

  ## Initialization
  step <- 1
  if (verbose) {
    print(paste("Step:", step))
    print(paste("|R1|=", 0))
  }
  thrMat <- NULL
  resJList <- list()
  converged <- FALSE

  ## joint FWER control through gammatification
  resJ <- getJointFWERThresholds(mat, refFamily="kFWER", ..., verbose=FALSE)
  resJList[[step]] <- resJ

  ##details<<At 'step 0', a one-parameter family \code{sLambda}
  ##ensuring JFWER control at level \code{alpha} is inferred using
  ##\code{getJointFWERThresholds}. This family is then kept fixed
  ##throughout the step-down process. What changes throughout the
  ##steps down is the over-estimation of \eqn{H_0}: starting with
  ##\eqn{H} at step 0, a first lower bound is derived from
  ##\code{sLambda(lambda)}, which allows re-calibration of
  ##\eqn{lambda} by applying \code{getJointFWERThresholds} to the
  ##subset of not-yet-rejected hypotheses, etc.

  ## One-parameter threshold family \eqn{s_lambda} in the BNR paper:
  sLambda <- resJ$sLambda
  ## Note that this family does not change during the step-down process!
  ## Only the value of \code{lambda} may change
  thr <- resJ$thr
  lambda <- resJ$lambda
  stopifnot(identical(thr, sLambda(lambda)))    ## sanity check
  thrMat <- cbind(thrMat, thr)

  ##details<<The default for \code{H0} is \code{NULL}, corresponding to the usual situation where the true null hypotheses are not known. The case where \code{H0} is not \code{NULL} implements an "Oracle" version of the step-down procedure, with only one step down.
  if (!is.null(H0)) {
      R1 <- setdiff(1:m, H0)
  } else {
      thr1 <- thr[1]   ## (1-)FWER threshold
      R1 <- which(stat>=thr1)
  }

  while (!converged) {
      step <- step+1
      if (verbose) {
          print(paste("Step:", step))
          print(paste("|R1|=", length(R1)))
      }

      if (length(R1)) {
          ## updated score matrix given R1
          if (verbose) {
              print(R1)
          }
          mat1 <- mat[-R1, ]
          sLambda1 <- function(alpha) sLambda(alpha)[-R1]

      } else {
          mat1 <- mat
          sLambda1 <- sLambda
      }

      ## joint FWER control through gammatification, *holding sLambda fixed*
      resJ <- getJointFWERThresholds(mat1, refFamily=sLambda1, ..., verbose=FALSE)
      resJList[[step]] <- resJ

      ## FWER threshold
      lambda <- resJ$lambda
      thr <- sLambda(lambda)
      thrMat <- cbind(thrMat, thr)

      if (!is.null(H0)) {
          R1new <- R1  ## force convergence
      } else {
          thr1 <- thr[1]   ## (1-)FWER threshold
          R1new <- which(stat>=thr1)
      }

      ## convergence reached?
      converged <- identical(R1new, R1)

      ## update R1
      R1 <- R1new
  }

  ##value<< List with elements:
  list(thr=thr, ##<< A numeric vector \code{thr} of length \eqn{m}, such that the estimated probability that there exists an index \eqn{k} between 1 and m such that the \eqn{k}-th maximum of the test statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}.
       thrMat=thrMat, ##<< A matrix of size \eqn{m} x the number of steps down, containing successive JFWER threshold estimates.
       steps=resJList ##<< A list of results of successive applications of 'getJointFWERThresholds'
       )
}, ex=function(){
  ## parameters
  m <- 1e3+1
  rho <- 0.2
  n <- 123
  pi0 <- 0.8
  B <- 1e3

  ##  set.seed(0xBEEF)
  sim <- simulateMein2006(m, rho, n, pi0, SNR=2)
  X <- sim$X
  y <- sim$y
  H0 <- which(sim$H==0)
  H1 <- which(sim$H==1)
  m0 <- length(H0)
  m1 <- length(H1)

  ## Test statistics
  w <- wilcoxStat(X, y, B=B)
  scoreMat <- w$stat0Mat
  stat <- w$stat

  if (FALSE) {
      ## show test statistics
      pch <- 20
      plot(stat, col=rep(c(1, 2), times=c(m0, m1)), main="Test statistics", pch=pch)
      legend("topleft", c("H0", "H1"), pch=pch, col=1:2)
  }
  alpha <- 0.1

  resSD <- stepDownControl(stat, scoreMat, refFamily="kFWER", alpha=alpha, verbose=TRUE)
  thrMat <- resSD$thrMat

  ## confidence envelopes
  nSteps <- ncol(thrMat)
  thr <- thrMat[, nSteps]
  o <- order(stat, decreasing=TRUE)
  statO <- stat[o]

  Vbar <- upperBoundFP(statO, thr)  ## default is flavor "Roquain2014"
  VbarM <- upperBoundFP(statO, thr, flavor="Mein2006")  ## faster for now
  identical(Vbar, VbarM)  ## Generally TRUE

  bounds <- apply(thrMat, 2, function(thr) upperBoundFP(statO, thr, flavor="Mein2006"))
  ## True V (number of false discoveries among first rejections)
  V <- cumsum(o %in% H0)

  ## comparison with "Oracle" step-down JFWER thresholds
  resO <- stepDownControl(stat, scoreMat, refFamily="kFWER", alpha=alpha, verbose=TRUE, H0=H0)  ## does this work?
  thrO <- resO$thr
  VbarO <- upperBoundFP(statO, thrO, flavor="Mein2006")

  ## comparison with "double Oracle" JFWER thresholds
  scoreMatOracle <- scoreMat[-H1, ]
  resO2 <- getJointFWERThresholds(scoreMatOracle, refFamily="kFWER", alpha=alpha)
  thrO2 <- c(resO2$thr, rep(-Inf, m1))
  VbarO2 <- upperBoundFP(statO, thrO2, flavor="Mein2006")

  cols1 <- seq.int(nSteps)
  ltys1 <- rep(1, nSteps)
  ttl <- paste("Bounds on #FP among rejected hypotheses",
               paste("m=", m, ", rho=", rho, ", alpha=", alpha, sep=""), sep="\n")
  xmax <- min(200, m)
  ymax <- bounds[xmax, 1]
  matplot(bounds, t='s', lty=1, ylab="V", col=cols1, main=ttl,
          xlim=c(1, xmax), ylim=c(0, ymax))

  cols2 <- c("purple", "pink", "orange")
  ltys2 <- rep(1, length(cols2))
  lines(V, col=cols2[1], t="s", lty=ltys2[1])
  lines(VbarO, col=cols2[2], t="s", lty=ltys2[2])
  lines(VbarO2, col=cols2[3], t="s", lty=ltys2[3])

  lgd <- c(paste("SD-JFWER(step=", 1:nSteps, ")", sep=""), "True V", "Oracle SD-JFWER", "Oracle JFWER")
  ltys <- c(ltys1, ltys2)
  cols <- c(cols1, cols2)
  legend("top", lgd, col=cols, lty=ltys)

  statSc <- stat/max(stat)*ymax  ## scaled test statistics (for display)
  colStat <- 1+sim$H
  points(statSc[o], col=colStat[o], pch=pch, cex=0.5)
  legend("left", c("H0", "H1"), pch=pch, col=1:2)
})

############################################################################
## HISTORY:
##
## 2016-01-07
## o BUG FIX: the one-parameter family should not be updated at each step down!
## o Example fixed accordingly.
## o Implemented the "Oracle" version of step-down control.
## 2014-05-09
## o Created from 'getJointFWERThresholds'.
############################################################################

