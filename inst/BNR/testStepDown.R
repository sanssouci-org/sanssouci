testStepDown <- function(m, dep, B, pi0, SNR, typeOfSNR, alphas, flavor=c("equi", "Mein2006", "Toeplitz", "equi-perm", "equi-flip"), kMaxs=m, trace=FALSE) {
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
    } else if (flavor=="equi-perm") {
        n <- 1e4 ## currently hardcoded (and small for speed reasons).
        p <- 0.5 ## currently hardcoded.
        rho <- dep
        sim <- simulateEquiByRandomization(m, rho, n, B, pi0, flavor="perm", p.value= FALSE, SNR=SNR, p=p, w=NULL)
    } else if (flavor=="equi-flip") {
        n <- 1e2 ## currently hardcoded (and small for speed reasons).
        rho <- dep
        sim <- simulateEquiByRandomization(m, rho, n, B, pi0, flavor="flip", p.value= FALSE, SNR=SNR, p=p, w=NULL)
    }
    X0 <- sim$X0
    x <- sim$x
    H <- sim$H

    ## Some data-driven rejection sets:
    ## subsample of a thresholding based rejection set
    pval <- 1-pnorm(x)  ## one-sided p-values
    wwBH5 <- userSelect(pval, 0.05, samplingFraction=1/2, method="BH")
    wwBH50 <- userSelect(pval, 0.5, samplingFraction=1/2, method="BH")
    ww0 <- userSelect(pval, 0.05, samplingFraction=1/2, method="none")
    rBH5 <- x[wwBH5]
    rBH50 <- x[wwBH50]
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
    H1 <- which(H==1)
    m0 <- length(H0)
    m1 <- m-m0
    nBH51 <- length(intersect(H1, wwBH5))
    nBH501 <- length(intersect(H1, wwBH50))
    n01 <- length(intersect(H1, ww0))
    
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
                rejkBH5 <- apply(resFam, 2, rejk, rBH5) ## subset of {p<BH(alpha)}
                rejkBH50 <- apply(resFam, 2, rejk, rBH50) ## subset of {p<BH(alpha)} 
                rejkP0 <- apply(resFam, 2, rejk, r0)  ## subset of {p<alpha}
                
                ## sanity checks
                stopifnot(all.equal(colSums(rejk0>=0), rej0))
                stopifnot(all.equal(colSums(rejk1>=0), rej1))
                
                v0 <- (m0-1-pmax(0, matrixStats::colMaxs(rejk0)))/m0  ## estimate of Vbar(H0)/m0: |H0|-1 - max_k |Rk \cup H0| - k =  min_k |Rk^c \cup H0| + (k-1) 
                s1 <- pmax(0, matrixStats::colMaxs(rejk1))/m1         ## estimate of Sbar(H1)/m1:          max_k |Rk \cap H1| - k
                ## s1b <- pmax(0, matrixStats::colMaxs(rejk1))/m1+1   ## estimate of Sbar(H1)/m1:          max_k |Rk \cap H1| - (k-1) ## =(by the BNR book)
                s01 <- pmax(0, matrixStats::colMaxs(rejk01))/m1       ## estimate of Sbar(H)/m1
                zBH5 <- pmax(0, matrixStats::colMaxs(rejkBH5))        ## estimate of Sbar(R_BH)
                zBH5 <- ifelse(nBH51==0, NA, zBH5/nBH51)
                zBH50 <- pmax(0, matrixStats::colMaxs(rejkBH50))        ## estimate of Sbar(R_BH)
                zBH50 <- ifelse(nBH501==0, NA, zBH50/nBH501)
                z0 <- pmax(0, matrixStats::colMaxs(rejkP0))           ## estimate of Sbar(R0)
                z0 <- ifelse(n01==0, NA, z0/n01)
                res <- rbind(JR=rej0>0, detPow1=rej1>0, detPow=rej01>0, v0, estPow1=s1, estPow=s01, powBH5=zBH5, powBH50=zBH50, pow0=z0)
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
    if (length(R) <= 1L) {  ## if length(R)==1L, then R is interpreted by 'sample' as the *number* of items to choose from :(
        ret <- R
    } else {
        probs <- rev(rank(pval[R]))/sum(1:length(R))
        ret <- try(sample(R, round(samplingFraction*length(R)), prob=probs))
        if (class(ret)=="try-error") {
            print(probs)
        }
    }
    ret
}
