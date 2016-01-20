m <- 1023
B <- 1e3

if (FALSE) {
    flavor <- c("independent", "equi-correlated", "3-factor model")[2]
    rho <- 0.2

    mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
} else {
    ## Toeplitz
    tcoefs <- toeplitz((1:m)^(-2))
    Sigma <- Matrix::Matrix(tcoefs, sparse = TRUE)
    mat <- simulateGaussianNullsFromSigma(m, B, Sigma)
}

str(mat)

alpha <- 0.2
thrMat <- NULL
lambdas <- NULL
probs <- NULL
cols <- NULL
ltys <- NULL

maxSteps <- 100
methods <- c("Simes", "kFWER")
kMaxs <- c(NA, 1, 2, 10, m)
for (mm in seq(along=methods)) {
    meth <- methods[mm]
    for (ss in seq(along=kMaxs)) {
        ms <- maxSteps
        kMax <- kMaxs[[ss]]
        if (is.na(kMax)) {
            ## just a trick to get the original thresholds
            kMax <- m
            ms <- 1  ## no dichotomy: original thresholds
        }
        res <- getJointFWERThresholds(mat, tau=meth, alpha, maxSteps=ms, kMax=kMax)
        thrMat <- cbind(thrMat, res$thr)
        lambdas <- c(lambdas, res$lambda)
        probs <- c(probs, res$prob)
        cols <- c(cols, ss)
        ltys <- c(ltys, mm)
    }
}
lgd <- paste(methods[ltys], "; k[max]=", kMaxs, "; lambda=", round(lambdas, 2), "; p=", round(probs, 2))
matplot(thrMat, t='l', log="x", col=cols, lty=ltys, cex=0.6)
matpoints(thrMat[1:10,], col=cols, cex=0.4, pch=1)
legend("bottomleft", as.expression(lgd), col=cols, lty=ltys)
