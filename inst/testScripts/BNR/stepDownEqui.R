library("sansSouci")


testStepDown <- function(m, rho, B, pi0, alpha) {
    sim <- simEqui(m, rho, B, pi0, SNR=SNR)
    X0 <- sim$X0
    x <- sim$x
    H <- sim$H
    
    ## truth
    H0 <- which(sim$H==0)
    H1 <- which(sim$H==1)
    m0 <- length(H0)
    m1 <- length(H1)
    cols <- c(1, 2)[1+H]
    
    # SD control
    resSD <- stepDownControl(x, X0, tau="kFWER", alpha=alpha, verbose=TRUE)
    thrMat <- resSD$thrMat
    
    # Final thresholds
    nSteps <- ncol(thrMat)
    thr <- thrMat[, nSteps]
    o <- order(x, decreasing=TRUE)
    xO <- x[o]
    
    # "Oracle" JFWER thresholds
    X0Oracle <- X0[-H1, ]
    resOracle <- getJointFWERThresholds(X0Oracle, tau="kFWER", alpha=alpha)
    thrO <- c(resOracle$thr, rep(-Inf, m1))
    
    # JFWER 
    nbBadK <- function(thr) {
        xNulls <- x[H0]
        nbFP <- sapply(thr, FUN=function(ss) sum(xNulls > ss))
        sum(nbFP >= 1:m)
    }
    
    allThr <- cbind(thrMat, thrO)
    badK <- apply(allThr, 2, nbBadK)
    return(badK)
}


# Parameters:
m <- 1e3
rho <- 0
pi0 <- 0.8
B <- 5e3
SNR <- 2
alpha <- 0.25

#set.seed(0xBEEF)
library("parallel")
nbSimu <- 8
ncores <- 4
res <- mclapply(1:nbSimu, FUN=function(bb) {
    print(bb)
    res <- testStepDown(m, rho, B, pi0, alpha)
    return(res)
}, mc.cores=ncores)

