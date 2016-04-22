library("sansSouci")
library("actuar")

Rcpp <- c(TRUE, FALSE)[1]

testStepDown <- function(m, rho, B, pi0, SNR, alpha, Rcpp=FALSE) {
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
    resSD <- stepDownControl(x, X0, refFamily="kFWER", alpha=alpha, Rcpp=Rcpp)
    thrMat <- resSD$thrMat
    nSteps <- ncol(thrMat)
    thrMat <- cbind("step0"=thrMat[, 1], "stepDown"=thrMat[, nSteps])

    ## comparison with *Oracle step-down* JFWER thresholds (only the step-down is Oracle)
    resOracleSD <- stepDownControl(x, X0, refFamily="kFWER", alpha=alpha, H0=H0, Rcpp=Rcpp)
    ##    thrOSD <- c(resOracleJ$thr, rep(-Inf, m1))
    thrOSD <- resOracleSD$thr

    ## *Oracle JFWER* thresholds (double Oracle!!!)
    X0Oracle <- X0[-H1, ]
    resOracleJ <- getJointFWERThresholds(X0Oracle, refFamily="kFWER", alpha=alpha, Rcpp=Rcpp)
    thrOJ <- c(resOracleJ$thr, rep(-Inf, m1))

    res <- cbind(x=x, truth=sim$H, thrMat, "Oracle"=thrOSD, "Oracle2"=thrOJ)
    return(res)
}


# Parameters:
m <- 1e3
rho <- 0
pi0 <- 0.9
B <- 1e4
#SNR <- 2
alpha <- 0.25
nbSimu <- 1e4

SNRs <- c(0, 1, 2, 3, 4, 5)
#SNRs <- "Pareto(2,2,1)"
#SNRs <- c("Pareto(0,2,1)", "Pareto(1,2,1)", "Pareto(2,2,1)"); library("actuar");
rhos <- c(0, 0.2, 0.4)
pi0s <- c(0.9, 0.99, 0.999)

library("parallel")
ncores <- 25

for (SNR in SNRs) {
    for (rho in rhos) {
        for (pi0 in pi0s) {


tags <- sprintf("m=%s,pi0=%s,rho=%s,SNR=%s,B=%s,alpha=%s,nbSimu=%s",
                m, pi0, rho, SNR, B, alpha, nbSimu)
if (Rcpp) {
    sname <- "stepDownEqui+power,Rcpp"
} else {
    sname <- "stepDownEqui+power"
}
            
filename <- sprintf("%s,%s.rds", sname, gsub("\\.", "_", tags))
path <- "resData"
path <- file.path(path, sname)
path <- R.utils::Arguments$getWritablePath(path)

print(tags)

set.seed(0xBEEF)
res <- mclapply(1:nbSimu, FUN=function(bb) {
    if (bb %% 100 == 0) {
        print(bb)
    }
    #print(system.time(res <- testStepDown(m, rho, B, pi0, alpha)))
    res <- testStepDown(m, rho, B, pi0, SNR, alpha)
    #str(res)
    return(res)
}, mc.cores=ncores)

pathname <- file.path(path, filename)
saveRDS(res, file=pathname)


#print(tags)
resMat <- Reduce(rbind, dat)
print(colMeans(resMat>0))

}}}
