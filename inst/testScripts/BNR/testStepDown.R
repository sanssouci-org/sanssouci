testStepDown <- function(m, rho, B, pi0, SNR, alpha, Rcpp=FALSE, maxTime=100) {
    pid <- Sys.getpid()
    t0 <- Sys.time()
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
    }
    sim <- simEqui(m, rho, B, pi0, SNR=SNR)
    X0 <- sim$X0
    x <- sim$x
    H <- sim$H

    t0f <-      format(t0, "%Y%m%d-%H:%M")
    filename <- sprintf("simTrace,%s,%s.rda", pid, t0f)
    save(x, X0, alpha=alpha, H=H, file=filename)

    ## m <- nrow(X0)
    ## B <- ncol(X0)
    H0 <- which(H==0)
    m0 <- length(H0)
    m1 <- m-m0

    resMat <- NULL
    for (refFam in c("kFWER", "Simes")) {
        res <- jointFWERControl(X0, refFamily=refFam, alpha=alpha, stat=x)
        thrMat <- res$stepsDown$thr
        nSteps <- ncol(thrMat)
        if (nSteps>2) print(nSteps)
        thr0 <- thrMat[, 1]
        thrSD <- res$thr
        stopifnot(identical(thrSD, thrMat[, nSteps])) ## sanity check

        ## comparison with *Oracle step-down* JFWER thresholds (only the step-down is Oracle)
        xx <- rep(Inf, length(x))
        xx[H0] <- -Inf
        resOracle <- jointFWERControl(X0, refFamily=refFam, alpha=alpha, stat=xx)
        thrO <- resOracle$thr
##        lambda <- resOracle$lambda
##        thr <- sk(lambda)

        ## *Oracle JFWER* thresholds (double Oracle!!!)
        X0Oracle <- X0[H0, ]
        resOracleJ <- jointFWERControl(X0Oracle, refFamily=refFam, alpha=alpha, maxStepsDown=1)  ## single step
        thrOJ <- c(resOracleJ$thr, rep(-Inf, m1))

        resFam <- cbind("0"=thr0, "SD"=thrSD, "Oracle"=thrO, "Oracle2"=thrOJ)
        colnames(resFam) <- paste(refFam, colnames(resFam), sep=".")
        resMat <- cbind(resMat, resFam)
    }
    file.remove(filename); ## delete trace file if the above code terminated
    if (FALSE) { ## too much disk space required!
        res <- cbind(x=x, truth=H, thrMat)
    } else {     ## summarize results
        x0 <- x[which(H==0)]
        x1 <- x[which(H==1)]
        rej0 <- apply(resMat, 2, rej, x0)
        rej1 <- apply(resMat, 2, rej, x1)
        res <- rbind(rej0, rej1)
    }
    return(res)
}

## to estimate JFWER and power
rej <- function(thr, x) {
    nover <- sapply(thr, FUN=function(ss) sum(x > ss))
    sum(nover >= 1:length(thr))
}

