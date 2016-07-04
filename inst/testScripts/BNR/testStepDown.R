testStepDown <- function(m, dep, B, pi0, SNR, typeOfSNR, alpha, flavor=c("equi", "Mein2006", "Toeplitz"), kMax=m, trace=FALSE) {
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

    if (trace) {
        t0f <-      format(Sys.time(), "%Y%m%d-%H:%M")
        filename <- sprintf("simTrace,%s,%s.rda", Sys.getpid(), t0f)
        save(x, X0, alpha=alpha, H=H, file=filename)
    }

    if (FALSE) {
        plot(x, col=1+H, main="Test statistics", pch=20)
    }

    ## m <- nrow(X0)
    ## B <- ncol(X0)
    H0 <- which(H==0)
    m0 <- length(H0)
    m1 <- m-m0

    resMat <- NULL
    for (refFam in c("kFWER", "Simes")) {
        ## Step-down control
        res <- jointFWERControl(X0, refFamily=refFam, alpha=alpha, stat=x, kMax=kMax, verbose=FALSE)
        thrMat <- res$stepsDown$thr
        nSteps <- ncol(thrMat)
        thr0 <- thrMat[, 1]
        thrSD <- res$thr
        stopifnot(identical(thrSD, thrMat[, nSteps])) ## sanity check

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

        resFam <- cbind("0"=thr0, "SD"=thrSD, "Oracle"=thrO, "Oracle2"=thrOJ)
        if (refFam=="Simes") {
            ## Simes control
            sk <- SimesThresholdFamily(m, kMax=kMax)
            resFam <- cbind(resFam, "alpha"=sk(alpha))
        }
        colnames(resFam) <- paste(refFam, colnames(resFam), sep=".")
        resMat <- cbind(resMat, resFam)
    }
    if (trace) {
        file.remove(filename); ## delete trace file if the above code terminated
    }
    if (FALSE) { ## too much disk space required!
        res <- cbind(x=x, truth=H, thrMat)
    } else {     ## summarize results
        x0 <- x[which(H==0)]
        x1 <- x[which(H==1)]
        rej0 <- apply(resMat, 2, rej, x0) ## nb of k s.t. |Rk \cup H0| > k-1 ('violators')
        rej1 <- apply(resMat, 2, rej, x1) ## nb of k s.t. |Rk \cup H1| > k-1 ('good catchers')

        rejk0 <- apply(resMat, 2, rejk, x0)
        rejk1 <- apply(resMat, 2, rejk, x1)

        ## sanity checks
        stopifnot(all.equal(colSums(rejk0>=0), rej0))
        stopifnot(all.equal(colSums(rejk1>=0), rej1))

        v0 <- (m0-1-pmax(0, matrixStats::colMaxs(rejk0)))/m0  ## estimate of Vbar(H0)/m0
        s1 <- pmax(0, matrixStats::colMaxs(rejk1))/m1         ## estimate of Sbar(H1)/m1

        res <- rbind(rej0, rej1, v0, s1)
    }
    return(res)
}

## to estimate JFWER and power
rej <- function(thr, x) {
    nover <- sapply(thr, FUN=function(ss) sum(x > ss))
    sum(nover >= 1:length(thr))
}

## to estimate JFWER and power
rejk <- function(thr, x) {
    nover <- sapply(thr, FUN=function(ss) sum(x > ss))
    nover-1:length(thr)
}
