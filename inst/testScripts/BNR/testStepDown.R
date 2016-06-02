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

    ## SD control
    resSD <- stepDownJointFWERControl(x, X0, refFamily="kFWER", alpha=alpha, Rcpp=Rcpp)
    thrMat <- resSD$thrMat
    nSteps <- ncol(thrMat)
    thrMat <- cbind("step0"=thrMat[, 1], "stepDown"=thrMat[, nSteps])

    ## comparison with *Oracle step-down* JFWER thresholds (only the step-down is Oracle)
    xx <- rep(Inf, length(x))
    xx[H0] <- -Inf
    resOracleSD <- stepDownJointFWERControl(xx, X0, refFamily="kFWER", alpha=alpha)
    thrOSD <- resOracleSD$thr

    ## *Oracle JFWER* thresholds (double Oracle!!!)
    X0Oracle <- X0[H0, ]
    resOracleJ <- jointFWERControl(X0Oracle, refFamily="kFWER", alpha=alpha, Rcpp=Rcpp)
    thrOJ <- c(resOracleJ$thr, rep(-Inf, m1))

    file.remove(filename); ## delete trace file if the above code worked!
    
    td <- Sys.time()-t0
    if (td>maxTime) {
        print(td)
        ##pryr::mem_used()
    }


    res <- cbind(x=x, truth=H, thrMat, "Oracle"=thrOSD, "Oracle2"=thrOJ)
    return(res)
}
