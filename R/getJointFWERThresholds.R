getJointFWERThresholds <- structure(function(
### get a family of thresholds that control the joint FWER
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test statistics under the null hypothesis. \describe{
### \item{m}{is the number of tested hypotheses}
### \item{B}{is the number of Monte-Carlo samples}}
    tau=c("Simes", "kFWER", "LR06"),
### A character value or a function that returns a vector of \code{m} thresholds (see \link{details}).
    alpha,
### Target joint-FWER level.
    kMax=nrow(mat),
### For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le k[max]}).
    maxSteps=100,
### Maximal number of steps in dichotomy.  If \code{maxSteps==1}, no dichotomy is performed
###
    tol=1e-4,
### Maximal tolerated distance between joint FWER level achieved by the
### result and the target level \eqn{\alpha}.
    verbose=FALSE
### If \code{TRUE}, print results of intermediate calculations.
### Defaults to \code{FALSE}.
    ) {
  m <- nrow(mat)
  B <- ncol(mat)
  if (alpha*B<1) {  ## sanity check
    stop("Please make sure that alpha*B is greater than 1")
  }
  matS <- apply(mat, 2, sort, decreasing=TRUE)
  ## truncate to [1,kMax]
  matS <- matS[1:kMax, , drop=FALSE]

  Q <- NULL
  if (mode(tau)=="character") {
###<<details{If \code{\tau} is a character value, it must be one of the following:\describe{
###   \item{Simes}{The classical family of thresholds introduced by Simes (1986): \eqn{\alpha*k/m}. This family yields joint FWER control if the test statistics are positively dependent (PRDS) under H0.}
###   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k, \eqn{(t_k)} controls the (marginal) k-FWER.}
### }
    tau <- match.arg(tau)
    if (tau=="kFWER") {
      etaMax <- min(2, 1/alpha)
      matSS <- apply(matS, 1, sort, decreasing=TRUE)
      matSS <- t(matSS)
      tau <- function(alpha) {
        idx <- floor(min(alpha,1)*(B-1))   ## NB: B-1 ensures *non-asymptotically* valid quantiles
        matSS[, idx]
      }
      Q <- matSS  ## Meinshausen's 'Q' matrix
    } else if (tau=="Simes") {
      etaMax <- min(20, 1/alpha)
      tau <- function(alpha) qnorm(1-min(alpha, 1)*(1:kMax)/m)
    }
  } else if (mode(tau)=="function") {
###<<details{If \code{\tau} is a function, it has to be of the form \eqn{\tau:\alpha \mapsto (\tau_k(\alpha))_{k=1 \dots m}}, such that \eqn{\tau(\alpha)} ensures marginal kFWER control, that is, \eqn{\forall k, P(k-inf(P_{i}) \leq \tau_k(alpha)) \leq \alpha}}
    etaMax <- 2
  }
  if (maxSteps==1) etaMax <- 2
  ## initialization
  etaMin <- 0;
  steps <- 0;
  prob <- 0;
  stopCrit <- FALSE;
  etas <- numeric(0);
  probs <- numeric(0);
  
  ## iteration
  while (!stopCrit && (steps<maxSteps)) {
    prob0 <- prob
    steps <- steps+1
    eta <- (etaMin+etaMax)/2  ## candidate value
    thr <- tau(alpha*eta)
    isAbove <- apply(matS, 2, ">", thr)
    ## details << A null hypothesis si rejected iff its test statistic is greater _or equal_ to a threshold value
    ## Hence the ">" and not ">=" in the definition of 'isAbove'
    dim(isAbove) <- c(length(thr), B); ## avoid coercion to vector if only one column (B==1)
    nAbove <- colSums(isAbove)
    prob <- mean(nAbove>0)
    if (verbose) {
      cat("Step:", steps, "\n")
      print(prob)
    }
    if (prob>alpha) {
      etaMax <- eta
    } else {
      etaMin <- eta      
    }
    etas <- c(etas, eta)
    probs <- c(probs, prob)
    
    admissible <- (alpha-prob>=0)
    converged <- admissible && ((alpha-prob)<=tol*alpha);
    ##    stuck <- admissible && (prob==prob0)  ## threshold family may change afterward, but only marginally so
    sameQuant <- (floor(etaMin*alpha*B)==floor(etaMax*alpha*B))
    stuck <- admissible && (prob==prob0 || sameQuant)
    stopCrit <- converged | stuck
  }
  if (steps==maxSteps) {
    reason <- "Maximal number of steps reached in dichotomy"
  } else if (stuck) {
    reason <- paste("Converged after", steps, "iterations without reaching target tolerance")
  } else {
    reason <- "Convergence reached"
  }
  if (!admissible) {
    warning("Could not achieve target JFWER control")
  }
  if (length(thr)==0) { ## may occur when 'idx' in 'tau' is 0...
    thr <- rep(Inf, kMax)
  }
  
  stopifnot(length(thr)==kMax)  ## sanity check
  idxs <- seq(from=kMax+1, to=m, length=m-kMax)
  thr[idxs] <- thr[kMax]
  stopifnot(length(thr)==m)  ## sanity check

  res <- list(thr=thr,
### A numeric vector \code{thr}, such that the estimated probability that
### there exists an index \eqn{k} between 1 and m such that the k-th maximum
### of the test statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}.
              prob=prob,
### Achieved proportion (should be between \eqn{alpha-tol} and \eqn{alpha}).
              probs=probs,
### The sequence of such proportions along the steps of the dichotomy.
              steps=steps,
### Number of dichotomy steps performed.
              lambda=eta,
### Correction factor (in [0,1]) from the original (kFWER) thresholds to JFWER thresholds.
              etas=etas,
### The sequence of such correction factors along the steps of the dichotomy.
              reason=reason,
### A character sequence, the reason for stopping.
              tau=tau,
### A function that returns a vector of \code{m} thresholds (see \link{details}).  It corresponds to the input argument \code{tau} if it was a function. Otherwise, it is calculated from the input matrix.
              Q=Q
### A \eqn{m} x \eqn{B} matrix of B realizations of ranked test statistics under H0
              )
}, ex=function(){
  m <- 1023
  B <- 1e3

  if (FALSE) {
    flavor <- c("independent", "equi-correlated", "3-factor model")[2]
    rho <- 0.2
    
    sim <- simulateGaussianNullsFromFactorModel(m, B=B, flavor=flavor, rho=rho)
  } else {
    ## Toeplitz
    tcoefs <- toeplitz((1:m)^(-2))
    Sigma <- Matrix(tcoefs, sparse = TRUE)
    sim <- simulateGaussianNullsFromSigma(m, B, Sigma)
  }
  
  mat <- sim$Y
  str(mat)
  image(sim$Sigma)

  alpha <- 0.2
  thrMat <- NULL
  lambdas <- NULL
  probs <- NULL
  cols <- NULL
  ltys <- NULL

  maxSteps <- 100
  methods <- c("Simes", "kFWER")
  kMaxs <- c(NA, 1, 2, 10, m)
  for (mm in seq(along=methods)) {
    meth <- methods[mm]
    for (ss in seq(along=kMaxs)) {
      ms <- maxSteps
      kMax <- kMaxs[[ss]]
      if (is.na(kMax)) {
        ## just a trick to get the original thresholds
        kMax <- m
        ms <- 1  ## no dichotomy: original thresholds
      }
      res <- getJointFWERThresholds(mat, tau=meth, alpha, maxSteps=ms, kMax=kMax)
      thrMat <- cbind(thrMat, res$thr)
      lambdas <- c(lambdas, res$lambda)
      probs <- c(probs, res$prob)
      cols <- c(cols, ss)
      ltys <- c(ltys, mm)
    }
  }
  lgd <- paste(methods[ltys], "; k[max]=", kMaxs, "; lambda=", round(lambdas, 2), "; p=", round(probs, 2))
  matplot(thrMat, t='l', log="x", col=cols, lty=ltys, cex=0.6)
  matpoints(thrMat[1:10,], col=cols, cex=0.4, pch=1)
  legend("bottomleft", as.expression(lgd), col=cols, lty=ltys)
})


############################################################################
## HISTORY:
## 2014-05-09
## o Added item 'tau' to return value.
## 2014-02-11
## o Now using '(B-1)' instead of 'B' in previous SPEEDUP.
## 2014-01-16
## o SPEEDUP: now pre-sorting 'matS' instead of using rowQuantiles within
## 'tau' for 'kFWER'.
## 2013-03-29
## o Created.
############################################################################

