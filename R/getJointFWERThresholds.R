getJointFWERThresholds <- structure(function(
### get a family of thresholds that control the joint FWER
    mat,
### A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
### statistics under the null hypothesis. \describe{ \item{m}{is the
### number of null hypotheses tested} \item{B}{is the number of
### Monte-Carlo samples}}
    refFamily=c("Simes", "kFWER"),
### A character value or a function that returns a vector of \code{m}
### thresholds (see details).
    alpha,
### Target joint FWER level.
    kMax=nrow(mat),
### For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
### k[max]}).
    flavor=c("dichotomy", "pivotalStat"),
### How should the threshold be calibrated? Defaults to
### 'dichotomy'. The other option ('pivotalStat') avoids dichotomy by
### directly calculating a pivotal statistic of which the target value
### of \eqn{lambda} is simply the quantile of order \eqn{alpha}.
    maxSteps=100,
### Maximal number of steps in dichotomy (hence only used when
### \code{flavor=='dichotomy'}).  If \code{maxSteps==1}, no dichotomy
### is performed.
    tol=1e-4,
### Maximal tolerated distance between joint FWER level achieved by the
### result and the target level \eqn{\alpha}.
    Rcpp=FALSE,
### If \code{TRUE}, some costly operations (sorting) are performed in C++.
    fullOutput=FALSE,
### Should marginal kFWERs and the equivalent of Meinshausen's matrix
### 'Q' be returned? Mainly used for troubleshooting.
    verbose=FALSE)
### If \code{TRUE}, print results of intermediate calculations.
### Defaults to \code{FALSE}.
    {
        flavor <- match.arg(flavor)
        if (verbose) {
            print(flavor)
        }
        m <- nrow(mat)
        B <- ncol(mat)
        if (alpha*B<1) {  ## sanity check
            stop("Please make sure that alpha*ncol(mat) is greater than 1!")
        }
        if (flavor=="pivotalStat") {
            warning("Flavor 'pivotalStat' not fully implemented yet!")
        }
        if (Rcpp) {
            kmaxH0 <- partialColSortDesc(mat, kMax);
        } else {
            ## k-max of the test statistics under H0:
            kmaxH0 <- apply(mat, 2, sort, decreasing=TRUE)
            ## truncate to [1,kMax]
            kmaxH0 <- kmaxH0[1:kMax, , drop=FALSE]
        }

        probk <- NULL;
        Q <- NULL;

        if (mode(refFamily)=="character") {
###<<details{If \code{\refFamily} is a character value, it must be one of the following:\describe{
###   \item{Simes}{The classical family of thresholds introduced by Simes (1986): \eqn{\alpha*k/m}. This family yields joint FWER control if the test statistics are positively dependent (PRDS) under H0.}
###   \item{kFWER}{A family \eqn{(t_k)} calibrated so that for each k, \eqn{(t_k)} controls the (marginal) k-FWER.}
### }
            refFamilyCh <- refFamily
            refFamily <- match.arg(refFamily)
            if (refFamily=="kFWER") {
                etaMax <- min(2, 1/alpha)
                if (Rcpp) {
                    Q <- rowSortDesc(kmaxH0)
                } else {
                    kmaxH0s <- apply(kmaxH0, 1, sort, decreasing=TRUE)
                    Q <- t(kmaxH0s) ## corresponds to matrix 'Q' in Meinshausen 2006
                }
                refFamily <- function(alpha) {
                    idx <- floor(min(alpha,1)*(B-1))   ## NB: B-1 ensures *non-asymptotically* valid quantiles
                    if (verbose) {
                        print(idx)
                    }
                    if (idx==0) {
                        thr <- rep(Inf, nrow(Q))
                    } else {
                        thr <- Q[, idx]
                    }
                    thr
                }
            } else if (refFamily=="Simes") {
                etaMax <- min(20, 1/alpha)
                refFamily <- function(alpha) qnorm(1-min(alpha, 1)*(1:kMax)/m)
            }
        } else if (mode(refFamily)=="function") {
###<<details{If \code{refFamily} is a function \eqn{tau}, it has to be of the form \eqn{tau:\alpha \mapsto (\tau_k(\alpha))_{k=1 \dots m}}, such that \eqn{\tau(\alpha)} ensures marginal kFWER control, that is, \eqn{\forall k, P(k-inf(P_{i}) \leq \tau_k(alpha)) \leq \alpha}}
            stopifnot(length(refFamily(alpha))==m)   ## sanity check
            refFamilyCh <- "user-defined"
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
                thr <- refFamily(alpha*eta)

                ## details << A null hypothesis si rejected iff its
                ## test statistic is greater _or equal_ to a threshold
                ## value. Hence the ">" and not ">=" in the below
                ## definition of 'isAbove'.
                if (Rcpp) {
                    prob <- coverage(thr, kmaxH0);
                    if (fullOutput) {
                        probk <- marginalKFWER(thr, kmaxH0);
                    }
                } else {
                    isAbove <- sweep(kmaxH0, 1, thr,  ">")
##                    browser()
                    dim(isAbove) <- c(length(thr), B); ## avoid coercion to vector if only one column (B==1)
                    nAbove <- colSums(isAbove)
                    prob <- mean(nAbove>0)
                    if (fullOutput) {
                        probk <- rowMeans(isAbove)
                    }
                }

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
            if (length(thr)==0) { ## may occur when 'idx' in 'refFamily' is 0...
                thr <- rep(Inf, kMax)
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

            stopifnot(refFamilyCh %in% c("Simes", "kFWER"))
            if (refFamilyCh=="Simes") {
            ## for Simes, s_k^{-1}(u) = (m/k)*(1-pnorm(u))
                pval <- 1-pnorm(kmaxH0)
                skInv <- sweep(pval, 1, m/1:m, "*")
                pivotalStat <- colMins(skInv)
            } else if (refFamilyCh=="kFWER") {
                rks <- rowRanks(kmaxH0) +1
                pivotalStat <- colMins(rks)/B
            }
            lambda <- quantile(pivotalStat, alpha, type=1)

            ## return values
            thr <- refFamily(lambda)
            if (Rcpp) {
                prob <- coverage(thr, kmaxH0);
                if (fullOutput) {
                    probk <- marginalKFWER(thr, kmaxH0);
                }
            } else {
                isAbove <- sweep(kmaxH0, 1, thr,  ">")
                dim(isAbove) <- c(length(thr), B); ## avoid coercion to vector if only one column (B==1)
                nAbove <- colSums(isAbove)
                prob <- mean(nAbove>0)
                if (fullOutput) {
                    probk <- rowMeans(isAbove)
                }
            }
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
        } ## if (flavor) { ...

        ##value<< List with elements:
        res <- list(thr=thr,  ##<< A numeric vector \code{thr}, such that the estimated probability that ##<< there exists an index \eqn{k} between 1 and m such that the k-th maximum ##<< of the test statistics of is greater than \eqn{thr[k]}, is less than \eqn{\alpha}.
                    prob=prob, ##<< Estimated probability that there exists an index \eqn{k} between 1 ##<< and m such that the k-th maximum of the test statistics of is ##<< greater than \eqn{thr[k]} (should be in \eqn{[\alpha-tol, alpha]}).
                    probs=probs, ##<< The sequence of such estimated probabilities along the steps of the dichotomy.
                    probk=probk, ##<< A vector of length \eqn{m} whose \eqn{k}th entry is the estimated probability that the k-th maximum of the test statistics of is greater than \eqn{thr[k]}. Only returned if argument \code{fullOutput} is set to \code{TRUE}.
                    steps=steps, ##<< Number of dichotomy steps performed.
                    lambda=lambda, ##<< JFWER threshold.
                    lambdas=lambdas, ##<< The sequence of candidate JFWER thresholds along the steps of dichotomy, or the JFWER threshold \eqn{lambda} if \code{flavor=="pivotalStat"}
                    reason=reason, ##<< A character sequence, the reason for stopping.
                    sLambda=refFamily, ##<< The function \eqn{lambda \mapsto sLambda}, such that \code{thr} is identical to \code{sLambda(lambda)}. If the input argument \code{refFamily} is a function, then \code{sLambda=refFamily}.
                    pivotalStat=pivotalStat, ##<< A numeric vector, the values of the pivotal statistic whose quantile of order \eqn{alpha} is \eqn{lambda}
                    ##end<<
                    Q=Q, ##<< Result of sorting the input score matrix by row and then by columns. Corresponds to matrix 'Q' in Meinshausen (2006).
                    kmaxH0=kmaxH0 ##<< k-max of the test statistics under H0 for \eqn{1 \leq k \leq kMax}.
#                    getCoverage=getCoverage  ##<< function that returns the 'coverage' of a threshold family in terms of JFWER.
                    )
### A \eqn{m} x \eqn{B} matrix of B realizations of ranked test statistics under H0
    }, ex=function(){
        m <- 1023
        B <- 1e3

        flavor <- c("independent", "equi-correlated", "3-factor model")[2]
        rho <- 0.2
        mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
        alpha <- 0.2

        res <- getJointFWERThresholds(mat, refFamily="kFWER", alpha)
        str(res)

        resP <- getJointFWERThresholds(mat, refFamily="kFWER", alpha, flavor="pivotalStat")
        str(resP)

    })


############################################################################
## HISTORY:
##
## 2016-01-28
## o Changed the input value from 'tau' to 'refFamily' for consistency
## with the notation of the BNR paper.
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

