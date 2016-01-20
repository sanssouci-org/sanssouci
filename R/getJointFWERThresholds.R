getJointFWERThresholds <- structure(function(
### get a family of thresholds that control the joint FWER
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
### statistics under the null hypothesis. \describe{ \item{m}{is the
### number of null hypotheses tested} \item{B}{is the number of
### Monte-Carlo samples}}
    tau=c("Simes", "kFWER"),
### A character value or a function that returns a vector of \code{m}
### thresholds (see \link{details}).
    alpha,
### Target joint-FWER level.
    kMax=nrow(mat),
### For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
### k[max]}).
    flavor=c("dichotomy", "pivotalStat"),
### How should the threshold be calibrated? Defaults to
### 'dichotomy'. The other option ('pivotalStat') avoids dichotomy by
### directly calculating a pivotal statistic of which the target value
### of \expr{lambda} is simply the quantile of order \expr{alpha}.
    maxSteps=100,
### Maximal number of steps in dichotomy (hence only used when
### \code{flavor=='dichotomy'}).  If \code{maxSteps==1}, no dichotomy
### is performed.
    tol=1e-4,
### Maximal tolerated distance between joint FWER level achieved by the
### result and the target level \eqn{\alpha}.
    verbose=FALSE)
### If \code{TRUE}, print results of intermediate calculations.
### Defaults to \code{FALSE}.
    {
        flavor <- match.arg(flavor)
        if (verbose) {
            print(flavor)
        }
        if (flavor=="pivotalStat" && B>1e3) {
            warning("The current implementation of flavor 'pivotalStat' is too slow for large values of 'B'. Forcing flavor='dichotomy'")
        }
        m <- nrow(mat)
        B <- ncol(mat)
        if (alpha*B<1) {  ## sanity check
            stop("Please make sure that alpha*ncol(mat) is greater than 1!")
        }
        ## k-max of the test statistics under H0:
        kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)
        ## truncate to [1,kMax]
        kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]

        Q <- NULL
        if (mode(tau)=="character") {
###<<details{If \code{\tau} is a character value, it must be one of the following:\describe{
###   \item{Simes}{The classical family of thresholds introduced by Simes (1986): \eqn{\alpha*k/m}. This family yields joint FWER control if the test statistics are positively dependent (PRDS) under H0.}
###   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k, \eqn{(t_k)} controls the (marginal) k-FWER.}
### }
            tauCh <- tau
            tau <- match.arg(tau)
            if (tau=="kFWER") {
                etaMax <- min(2, 1/alpha)
                kmaxH0s <- apply(kmaxH0, 1, sort, decreasing=TRUE)
                Q <- t(kmaxH0s) ## corresponds to matrix 'Q' in Meinshausen 2006
                tau <- function(alpha) {
                    idx <- floor(min(alpha,1)*(B-1))   ## NB: B-1 ensures *non-asymptotically* valid quantiles
                    Q[, idx]
                }
            } else if (tau=="Simes") {
                etaMax <- min(20, 1/alpha)
                tau <- function(alpha) qnorm(1-min(alpha, 1)*(1:kMax)/m)
            }
        } else if (mode(tau)=="function") {
###<<details{If \code{\tau} is a function, it has to be of the form \eqn{\tau:\alpha \mapsto (\tau_k(\alpha))_{k=1 \dots m}}, such that \eqn{\tau(\alpha)} ensures marginal kFWER control, that is, \eqn{\forall k, P(k-inf(P_{i}) \leq \tau_k(alpha)) \leq \alpha}}
            stopifnot(tau(alpha)==m)   ## sanity check
            tauCh <- "user-defined"
            etaMax <- 2
        }

        if (flavor=="dichotomy") {

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
                isAbove <- sweep(kmaxH0, 1, thr,  ">")
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
            } else {
                probk <- rowMeans(isAbove)
            }
            pivotalStat <- NA_real_
            lambdas <- alpha*etas
            lambda <- alpha*eta

            stopifnot(length(thr)==kMax)  ## sanity check
            idxs <- seq(from=kMax+1, to=m, length=m-kMax)
            thr[idxs] <- thr[kMax]
            stopifnot(length(thr)==m)  ## sanity check
        } else {            ## flavor=="pivotalStat"
            ###<<details{If \code{flavor=="pivotalStat"}, then the
            ###function returns an element \code{pivotalStat}, whose
            ###quantile of order \code{alpha} is the target
            ###\code{lambda(\alpha)}. A major advantage of this option
            ###is that \lambda(\alpha) may then be calculated for any
            ###\code{\alpha} at no additional comutational cost. This
            ###flavor is slightly faster than the default flavor for
            ###Simes thresholds, but it is substantially slower (with
            ###the current implementation) for kFWER thresholds.}

            stopifnot(tauCh %in% c("Simes", "kFWER"))
            if (tauCh=="Simes") {
            ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
                pval <- 1-pnorm(kmaxH0)
                skInv <- sweep(pval, 1, m/1:m, "*")
                mins <- colMins(skInv)
            } else if (tauCh=="kFWER") {
                getMin <- function(bb, Q) {
                    kb <- kmaxH0[,bb]
                    sw <- sweep(Q, 1, kb, "<")
                    Fkhat <- rowSums(sw)
                    return(min(Fkhat))
                }
                mins <- sapply(1:B, getMin, Q)/B   ## B+1???
            }
            lambda <- quantile(mins, alpha, type=1)

            ## return values
            thr <- tau(lambda)
            isAbove <- sweep(kmaxH0, 1, thr,  ">")
            dim(isAbove) <- c(length(thr), B); ## avoid coercion to vector if only one column (B==1)
            nAbove <- colSums(isAbove)
            probk <- rowMeans(isAbove)
            prob <- mean(nAbove>0)
            if (verbose) {
                print(prob)
            }
            if (prob>alpha) {
                stop("'prob' should be less than 'alpha'!")
            }
            probs <- prob
            steps <- 0L
            lambdas <- lambda
            reason <- NA_character_
            pivotalStat <- mins
        } ## if (flavor) { ...

        ##value<< List with elements:
        res <- list(thr=thr,  ##<< A numeric vector \code{thr}, such that the estimated probability that ##<< there exists an index \eqn{k} between 1 and m such that the k-th maximum ##<< of the test statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}.
                    prob=prob, ##<< Estimated probability that there exists an index \eqn{k} between 1 ##<< and m such that the k-th maximum of the test statistics of is ##<< greater than \eqn{thr[k]} (should be in \eqn{[\alpha-tol, alpha]}).
                    probs=probs, ##<< The sequence of such estimated probabilities along the steps of the dichotomy.
                    probk=probk, ##<< A vector of length \eqn{m} whose \eqn{k}th entry is the estimated probability that the k-th maximum of the test statistics of is greater than \eqn{thr[k]}.
                    steps=steps, ##<< Number of dichotomy steps performed.
                    lambda=lambda, ##<< JFWER threshold.
                    lambdas=lambdas, ##<< The sequence of candidate JFWER thresholds along the steps of dichotomy, or the JFWER threshold \eqn{lambda} if \code{flavor=="pivotalStat"}
                    reason=reason, ##<< A character sequence, the reason for stopping.
                    tau=tau, ##<< A function that returns a vector of \code{m} thresholds (see \link{details}).  It corresponds to the input argument \code{tau} if it was a function. Otherwise, it is calculated from the input matrix.
                    Q=Q, ##<< Result of sorting the input score matrix by row and then by columns. Corresponds to matrix 'Q' in Meinshausen (2006).
                    sLambda=tau, ##<< The function \eqn{lambda \mapsto s(lambda)}, such that \code{thr} is identical to \code{sLambda(lambda)}.
                    pivotalStat=pivotalStat ##<< A numeric vector, the values of the pivotal statistic whose quantile of order \eqn{alpha} is \eqn{lambda}
                    ##end<<
                    )
### A \eqn{m} x \eqn{B} matrix of B realizations of ranked test statistics under H0
    }, ex=function(){
        m <- 1023
        B <- 1e3

        flavor <- c("independent", "equi-correlated", "3-factor model")[2]
        rho <- 0.2
        mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
        alpha <- 0.2

        res <- getJointFWERThresholds(mat, tau="kFWER", alpha)
        str(res)

        resP <- getJointFWERThresholds(mat, tau="kFWER", alpha, flavor="pivotalStat")
        str(resP)

    })


############################################################################
## HISTORY:
##
## 2016-01-07
## o Changed the return value from 'eta' to 'lambda' for consistency
## with the notation of the BNR paper.
## o Added 'sLambda' to the return values of the function.
## 2015-12-10
## o Added argument 'flavor' to make use of the pivotal statistic
## identified by Etienne (and thus avoid dichotomy).
## 2014-05-09
## o Added item 'tau' to return value.
## 2014-02-11
## o Now using '(B-1)' instead of 'B' in previous SPEEDUP.
## 2014-01-16
## o SPEEDUP: now pre-sorting 'kmaxH0' instead of using rowQuantiles within
## 'tau' for 'kFWER'.
## 2013-03-29
## o Created.
############################################################################

