stepDownControl <- structure(function(
### get a family of thresholds that provide step down control of the joint FWER
    stat,
### A vector of \eqn{m} test statistics
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test statistics under the null hypothesis. \describe{
### \item{m}{is the number of tested hypotheses}
### \item{B}{is the number of Monte-Carlo samples}}
    tau=c("Simes", "kFWER", "LR06"),
### A character value or a function that returns a vector of \code{m} thresholds (see \link{details}).
    ...,
### Further arguments to be passed to getJointFWERControl.
    verbose=FALSE
### If \code{TRUE}, print results of intermediate calculations.
### Defaults to \code{FALSE}.
    ) {
  m <- nrow(mat)
  stopifnot(length(stat)==m)
  stopifnot(tau=="kFWER")  ## other flavor not implemented yet (?)
  
  ## Generic 'tau' function
  tauQR <- function(alpha, Q, r) {
    m <- nrow(Q)
    B <- ncol(Q)
    idx <- floor(min(alpha,1)*(B-1))   ## NB: B-1 ensures *non-asymptotically* valid quantiles
    Q[1:(m-r), idx]  ## A vector \eqn{\tau(\alpha)} which ensures marginal kFWER control, that is, \eqn{\forall k, P(k-inf(P_{i}) \leq \tau_k(alpha)) \leq \alpha}}
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
  resJ <- getJointFWERThresholds(mat, tau="kFWER", ..., verbose=FALSE)
  resJList[[step]] <- resJ

  thr <- resJ$thr
  thrMat <- cbind(thrMat, thr)

  thr1 <- thr[1]  ## FWER threshold
  R1 <- which(stat>=thr1)
  r1 <- length(R1)
  
  ## updated threshold function after removing |R1| smallest items
  Q <- resJ$Q
  stopifnot(!is.null(Q)) ## sanity check: 'Q' is not NULL (expecting a matrix...)
  tau <- function(alpha) tauQR(alpha, Q, r1)

  while (!converged) {
    step <- step+1
    if (verbose) {
      print(paste("Step:", step))
      print(paste("|R1|=", r1))
    }
    
    if (length(R1)) {
      ## updated score matrix given R1
      mat1 <- mat[-R1, ]
    } else {
      mat1 <- mat
    }
    
    ## joint FWER control through gammatification
    resJ <- getJointFWERThresholds(mat1, tau=tau, ..., verbose=FALSE)
    resJList[[step]] <- resJ
    
    ## FWER threshold
    thr <- resJ$thr
    R1new <- which(stat>=thr[1])

    ## fill the vector of joint FWER thresholds with -Inf (if k>|H_0|, R_k=H works)
    thr <- c(thr, rep(-Inf, length(R1)))
    stopifnot(length(thr)==m)               ## sanity check
    thrMat <- cbind(thrMat, thr)

    ## convergence reached?
    converged1 <- (length(R1new)==r1)
    converged <- identical(R1new, R1)
    stopifnot(converged==converged1)

    ## update R1 and friends
    R1 <- R1new
    r1 <- length(R1)
    
    ## updated threshold function after removing |R1| smallest items
    tau <- function(alpha) tauQR(alpha, Q, r1)
  }

  list(thr=thr,
### A numeric vector \code{thr} of length \eqn{m}, such that the estimated probability that
### there exists an index \eqn{k} between 1 and m such that the k-th maximum
### of the test statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}.
       thrMat=thrMat,
### A matrix of size \eqn{m} x the number of steps down, containing
### successive JFWER threshold estimates.
       steps=resJList
### A list of results of successive applications of 'getJointFWERThresholds'
       )
}, ex=function(){
  ## parameters
  m <- 1e3
  rho <- 0.2
  n <- 123
  pi0 <- 0.8
  B <- 1e4
  
  set.seed(0xBEEF)
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

  ## show test statistics
  pch <- 20
  plot(stat, col=rep(c(1, 2), times=c(m0, m1)), main="Test statistics", pch=pch)
  legend("topleft", c("H0", "H1"), pch=pch, col=1:2)

  alpha <- 0.1

  resSD <- stepDownControl(stat, scoreMat, tau="kFWER", alpha=alpha, verbose=TRUE)
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

  ## comparison with "Oracle" JFWER thresholds
  scoreMatOracle <- scoreMat[-H1, ]
  resOracle <- getJointFWERThresholds(scoreMatOracle, tau="kFWER", alpha=alpha)
  thrO <- c(resOracle$thr, rep(-Inf, m1))
  VbarO <- upperBoundFP(statO, thrO, flavor="Mein2006")
  
  cols1 <- seq.int(nSteps)
  ltys1 <- rep(1, nSteps)
  ttl <- paste("Bounds on #FP among rejected hypotheses",
               paste("m=", m, ", rho=", rho, ", alpha=", alpha, sep=""), sep="\n")
  xmax <- min(200, m)
  ymax <- bounds[xmax, 1]
  matplot(bounds, t='s', lty=1, ylab="V", col=cols1, main=ttl,
          xlim=c(1, xmax), ylim=c(0, ymax))

  cols2 <- c("purple", "pink")
  ltys2 <- rep(1, 2)
  lines(V, col=cols2[1], t="s", lty=ltys2[1])
  lines(VbarO, col=cols2[2], t="s", lty=ltys2[2])
  
  lgd <- c(paste("SD-JFWER(step=", 1:nSteps, ")", sep=""), "True V", "Oracle JFWER")
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
## 2014-12-10
## o 'tauQR' depends on 'R' only through its length.
## 
## 2014-05-09
## o Created from 'getJointFWERThresholds'.
############################################################################

