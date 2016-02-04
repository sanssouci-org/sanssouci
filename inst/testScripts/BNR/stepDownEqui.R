library("sansSouci")
library("actuar")

testStepDown <- function(m, rho, B, pi0, SNR, alpha) {
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

    ## truth
    H0 <- which(sim$H==0)
    H1 <- which(sim$H==1)
    m0 <- length(H0)
    m1 <- length(H1)
#    cols <- c(1, 2)[1+H]

    # SD control
    resSD <- stepDownControl(x, X0, refFamily="kFWER", alpha=alpha)
    thrMat <- resSD$thrMat

    # Final thresholds
    nSteps <- ncol(thrMat)
    thr <- thrMat[, nSteps]
    o <- order(x, decreasing=TRUE)

    ## comparison with *Oracle step-down* JFWER thresholds (only the step-down is Oracle)
    resOracleSD <- stepDownControl(x, X0, refFamily="kFWER", alpha=alpha, H0=H0)
    ##    thrOSD <- c(resOracleJ$thr, rep(-Inf, m1))
    thrOSD <- resOracleSD$thr

    ## *Oracle JFWER* thresholds (double Oracle!!!)
    X0Oracle <- X0[-H1, ]
    resOracleJ <- getJointFWERThresholds(X0Oracle, refFamily="kFWER", alpha=alpha)
    thrOJ <- c(resOracleJ$thr, rep(-Inf, m1))

    # JFWER
    nbBadK <- function(thr) {
        xNulls <- x[H0]
        nbFP <- sapply(thr, FUN=function(ss) sum(xNulls > ss))
        sum(nbFP >= 1:m)
    }

    allThr <- cbind(thrMat, thrOSD, thrOJ)
    badK <- apply(allThr, 2, nbBadK)
    return(badK)
}


# Parameters:
m <- 1e3
rho <- 0
pi0 <- 0.9
B <- 5e3
#SNR <- 2
alpha <- 0.25
nbSimu <- 1e2

SNRs <- c(0, 1, 2, 3, 4, 5)
#SNRs <- "Pareto(2,2,1)"
rhos <- c(0, 0.2, 0.4)
pi0s <- c(0.9, 0.99, 0.999)

library("parallel")
ncores <- 4

for (SNR in SNRs) {
    for (rho in rhos) {
        for (pi0 in pi0s) {


tags <- sprintf("m=%s,pi0=%s,rho=%s,SNR=%s,B=%s,alpha=%s,nbSimu=%s",
                m, pi0, rho, SNR, B, alpha, nbSimu)
sname <- "stepDownEqui"
filename <- sprintf("%s,%s.rds", sname, gsub("\\.", "_", tags))
path <- "resData"
path <- file.path(path, sname)
path <- R.utils::Arguments$getWritablePath(path)

print(tags)

#set.seed(0xBEEF)
res <- mclapply(1:nbSimu, FUN=function(bb) {
    if (bb %% 100 == 0) {
        print(bb)
    }
    #print(system.time(res <- testStepDown(m, rho, B, pi0, alpha)))
    res <- testStepDown(m, rho, B, pi0, SNR, alpha)
    return(res)
}, mc.cores=ncores)

pathname <- file.path(path, filename)
saveRDS(res, file=pathname)


#print(tags)
resMat <- sapply(res, FUN=function(x) x[c(1, length(x)-c(1, 0))])
rowMeans(resMat>0)

}}}
