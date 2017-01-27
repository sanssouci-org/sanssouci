testStepDown <- function(m, dep, B, pi0, SNR, typeOfSNR, alphas, flavor=c("equi", "Mein2006", "Toeplitz"), kMaxs=m, trace=FALSE) {
    flavor <- match.arg(flavor)
    if (is.character(SNR)) {
        pattern <- "Pareto\\(([0-9]),([0-9]),([0-9]))"
        xmin <- as.numeric(gsub(pattern, "\\1", SNR))
        shape <- as.numeric(gsub(pattern, "\\2", SNR))
        scale <- as.numeric(gsub(pattern, "\\3", SNR))

        if (is.na(xmin) || is.na(shape) || is.na(scale)) {
            stop("'SNR' parameter ill-specified in 'testStepDown'")
        }

        m1 <- round(m*(1-pi0))
        SNR <- xmin + rpareto(m1, shape, scale)
    } else if (typeOfSNR=="linear") {
        m1 <- round(m*(1-pi0))
        SNR <- seq(from=0, to=SNR, length=m1)
    }


    if (flavor=="equi") {
        rho <- dep
        parName <- "rho"
        sim <- simulateEqui(m, rho, B, pi0, SNR=SNR)
    } else if (flavor=="Mein2006") {
        n <- 1e3 ## currently hardcoded.
        rho <- dep
        sim <- simulateMein2006(m, rho, n, B, pi0, SNR=SNR)
    } else if (flavor=="Toeplitz") {
        ## Toeplitz, long range
        pow <- dep
        parName <- "pow"
        sim <- simulateToeplitz(m, pow, B, pi0, SNR=SNR)
    }
    X0 <- sim$X0
    x <- sim$x
    H <- sim$H

    ## Some data-driven rejection sets:
    ## subsample of a thresholding based rejection set
    pval <- 1-pnorm(x)  ## one-sided p-values
    wwBH <- userSelect(pval, 0.05, samplingFraction=1/2, method="BH")
    ww0 <- userSelect(pval, 0.05, samplingFraction=1/2, method="none")
    rBH <- x[wwBH]
    r0 <- x[ww0]
    
    if (trace) {
        t0f <-      format(Sys.time(), "%Y%m%d-%H:%M")
        filename <- sprintf("simTrace,%s,%s.rda", Sys.getpid(), t0f)
        save(x, X0, H=H, file=filename)
    }

    if (FALSE) {
        plot(x, col=1+H, main="Test statistics", pch=20)
    }

    ## m <- nrow(X0)
    ## B <- ncol(X0)
    H0 <- which(H==0)
    m0 <- length(H0)
    m1 <- m-m0

    resList <- NULL
    for (alpha in alphas) {
        ## atag <- sprintf("alpha=%s", alpha)
        atag <- as.character(alpha)
        for (kMax in kMaxs) {
            ## ktag <- sprintf("kMax=%s", kMax)
            ktag <- as.character(kMax)
            for (refFam in c("kFWER", "Simes")) {
                ## ftag <- sprintf("family=%s", refFam)
                ftag <- refFam
                ## Step-down control
                res <- jointFWERControl(X0, refFamily=refFam, alpha=alpha, stat=x, kMax=kMax, verbose=FALSE)
                thrMat <- res$stepsDown$thr
                nSteps <- ncol(thrMat)
                thr0 <- thrMat[, 1]
                thrSD <- res$thr
#                stopifnot(identical(thrSD, thrMat[, nSteps])) ## sanity check

                ## comparison with *Oracle step-down* JFWER thresholds (only the step-down is Oracle)
                xx <- rep(Inf, length(x))
                xx[H0] <- -Inf
                resOracle <- jointFWERControl(X0, refFamily=refFam, alpha=alpha, stat=xx, kMax=kMax, verbose=FALSE)
                thrO <- resOracle$thr

                ## *Oracle JFWER* thresholds (double Oracle!!!)
                X0Oracle <- X0[H0, ]
                resOracleJ <- jointFWERControl(X0Oracle, refFamily=refFam, alpha=alpha,
                                               maxStepsDown=0, kMax=kMax, verbose=FALSE)  ## single step
                thrOJ <- resOracleJ$thr
                if (m0<kMax) {
                    thrOJ <- c(thrOJ, rep(-Inf, kMax-m0))
                }

                resFam <- cbind("Single Step"=thr0, "Step down"=thrSD, "Oracle"=thrO, "Oracle2"=thrOJ)
                if (refFam=="Simes") {
                    ## Simes control
                    sk <- SimesThresholdFamily(m, kMax=kMax)
                    resFam <- cbind(resFam, "unadjusted"=sk(alpha))
                } else if (refFam=="kFWER") {
                    ## marginal kFWER control
                    sk <- kFWERThresholdFamily(X0, kMax=kMax)
                    resFam <- cbind(resFam, "unadjusted"=sk(alpha))
                }

                if (trace) {
                    file.remove(filename); ## delete trace file if the above code terminated
                }

                ## summarize results (o/w would require too much disk space)
                x0 <- x[which(H==0)]
                x1 <- x[which(H==1)]

                rej0 <- apply(resFam, 2, rej, x0) ## nb of k / |Rk \cup H0| > k-1 ('violators')
                rej1 <- apply(resFam, 2, rej, x1) ## nb of k / |Rk \cup H1| > k-1 ('good catchers')
                rej01 <- apply(resFam, 2, rej, x) ## nb of k / |Rk \cup H|  > k-1 ('detections')
                
                rejk0 <- apply(resFam, 2, rejk, x0)
                rejk1 <- apply(resFam, 2, rejk, x1)
                rejk01 <- apply(resFam, 2, rejk, x)
                rejkBH <- apply(resFam, 2, rejk, rBH) ## subset of {p<BH(alpha)} 
                rejkP0 <- apply(resFam, 2, rejk, r0)  ## subset of {p<alpha}
                
                ## sanity checks
                stopifnot(all.equal(colSums(rejk0>=0), rej0))
                stopifnot(all.equal(colSums(rejk1>=0), rej1))
                
                v0 <- (m0-1-pmax(0, matrixStats::colMaxs(rejk0)))/m0  ## estimate of Vbar(H0)/m0: |H0|-1 - max_k |Rk \cup H0| - k =  min_k |Rk^c \cup H0| + (k-1) 
                s1 <- pmax(0, matrixStats::colMaxs(rejk1))/m1         ## estimate of Sbar(H1)/m1:          max_k |Rk \cap H1| - k
                ## s1b <- pmax(0, matrixStats::colMaxs(rejk1))/m1+1   ## estimate of Sbar(H1)/m1:          max_k |Rk \cap H1| - (k-1) ## =(by the BNR book)
                s01 <- pmax(0, matrixStats::colMaxs(rejk01))/m1       ## estimate of Sbar(H)/m1
                zBH <- pmax(0, matrixStats::colMaxs(rejkBH))/m1       ## estimate of Sbar(R_BH)/m1
                z0 <- pmax(0, matrixStats::colMaxs(rejkP0))/m1        ## estimate of Sbar(R0)/m1                
                res <- rbind(JR=rej0>0, detPow1=rej1>0, detPow=rej01>0, v0, estPow1=s1, estPow=s01, powBH=zBH, pow0=z0)
                resList[[atag]][[ktag]][[ftag]] <- res
                rm(resFam, res);
            }
        }
    }
    return(resList)
}

## to estimate JFWER and power
rej <- function(thr, x) {
    nover <- sapply(thr, FUN=function(ss) sum(x > ss))  ## nover[k] is the number of items in x not rejected by R_k
    sum(nover >= 1:length(thr))                         ## number of k:s for which nover[k]>=k. If x=H0, this is the number of violators
}

## to estimate JFWER and power
rejk <- function(thr, x) {
    nover <- sapply(thr, FUN=function(ss) sum(x > ss))  ## nover[k] is the number of items in x not rejected by R_k
    nover-1:length(thr)                                 ## nover[k]-k. If x=H0, positive items correspond to violators
}

## subsample of thresholding-based rejection set
userSelect <- function(pval, alpha, samplingFraction=0.5, method="BH") {
    R <- which(p.adjust(pval, method=method)<alpha)  ## uninteresting because too adaptative ?
    if (length(R)==1L) {
        ret <- R
    } else {
        probs <- rev(rank(pval[R]))/sum(1:length(R))
        ret <- sample(R, round(samplingFraction*length(R)), prob=probs)
    }
    ret
}