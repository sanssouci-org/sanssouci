library("sansSouci")
library("parallel")

## model parameters
m <- 1e3
kMax <- 100
pi0 <- 0.9

## for (kMax in  c(m, 100)) {
##     for (pi0 in  c(0.8, 0.9, 0.95, 0.99)) {

B <- 5e3+1

## JFWER setup
maxSteps <- 100

flavor <- "equi-correlated"
rhos <- seq(from=0, to=0.9, by=0.1)

alpha <- 0.1

m0 <- round(m*pi0)
m1 <- m-m0
H0 <- sort(sample(m, m0))

rhos <- seq(from=0, to=0.5, by=0.1)

resK <- mclapply(1:length(rhos), FUN=function(rr) {
    rho <- rhos[rr]
    print(rho)

    mat <- simulateGaussianNullsFromFactorModel(m, n=B, flavor=flavor, rho=rho)
    mat0 <- mat[H0, ]     ## to mimic oracle thresholds for step-down

    ## Simes
    res <- getJointFWERThresholds(mat, tau="Simes", alpha, maxSteps=maxSteps, kMax=kMax)
    pk <- c(res$probk, rep(NA, m-length(res$probk)))
    mMax <- min(kMax, m0)
    sLambda0 <- function(alpha) res$sLambda(alpha)[1:mMax]
    res0 <- getJointFWERThresholds(mat0, tau=sLambda0, alpha, maxSteps=maxSteps, kMax=mMax)
    pk0 <- c(res0$probk, rep(NA, m-length(res0$probk)))
    pSimes <- cbind(pk, pk0)
    colnames(pSimes) <- c("Simes", "Simes H0")

    ## kFWER
    res <- getJointFWERThresholds(mat, tau="kFWER", alpha, maxSteps=maxSteps, kMax=kMax)
    pk <- c(res$probk, rep(NA, m-length(res$probk)))
    sLambda0 <- function(alpha) res$sLambda(alpha)[1:mMax]
    res0 <- getJointFWERThresholds(mat0, tau=sLambda0, alpha, maxSteps=maxSteps, kMax=mMax)
    pk0 <- c(res0$probk, rep(NA, m-length(res0$probk)))
    pKFWER <- cbind(pk, pk0)
    colnames(pKFWER) <- c("kFWER", "kFWER H0")

    return(cbind(pSimes, pKFWER))
}, mc.cores=3)

resSimes <- sapply(resK, FUN=function(x) x[, "Simes"])
resSimesH0 <- sapply(resK, FUN=function(x) x[, "Simes H0"])

resKFWER <- sapply(resK, FUN=function(x) x[, "kFWER"])
resKFWERH0 <- sapply(resK, FUN=function(x) x[, "kFWER H0"])

colnames(resSimes) <- rhos
colnames(resKFWER) <- rhos

ww <- match(c(0, 0.1, 0.2, 0.5), rhos)
cols <- seq(along=ww)

tags <- sprintf("m=%s,m0=%s,kMax=%s", m, m0, kMax)
figPath <- "fig"
dir.create(figPath, showWarnings=FALSE)

figName <- "unbalanced"
filename <- sprintf("%s,%s.pdf", figName, tags)
pathname <- file.path(figPath, filename)
pdf(pathname)
par(lwd=1.5)
matplot(resSimes[, ww], t='l', lty=1, col=cols, ylab=expression(P(p[(k)] < tau[k])), xlab="k", xlim=c(1, kMax))
matlines(resKFWER[, ww], t='l', lty=2, col=cols)
matlines(resKFWERH0[, ww], t='l', lty=4, col=cols)
legend("topright", col=cols, as.character(rhos[ww]), title=expression(rho), lty=1)
legend("topleft", col=1, lty=c(1:2, 4), c("Simes", "kFWER", "kFWER-H0"), horiz=TRUE)
dev.off()


figName <- "unbalanced+0"
filename <- sprintf("%s,%s.pdf", figName, tags)
pathname <- file.path(figPath, filename)
pdf(pathname)
par(lwd=1.5)
matplot(resSimes[, ww], t='l', lty=1, col=cols, ylab=expression(P(p[(k)] < tau[k])), xlab="k", xlim=c(1, kMax))
matlines(resKFWER[, ww], t='l', lty=2, col=cols)
matlines(resKFWERH0[, ww], t='l', lty=4, col=cols)
legend("topright", col=cols, as.character(rhos[ww]), title=expression(rho), lty=1)
legend("topleft", col=1, lty=c(1:2, 4), c("Simes", "kFWER", "kFWER-H0"), horiz=TRUE)
dev.off()

figName <- "unbalanced,zoom"
filename <- sprintf("%s,%s.pdf", figName, tags)
pathname <- file.path(figPath, filename)
pdf(pathname)
matplot(resSimes[, ww], t='l', lty=1, col=cols, ylab=expression(P(p[(k)] < tau[k])), xlim=c(1,20), xlab="k")
matlines(resKFWER[, ww], t='l', lty=2, col=cols)
matlines(resKFWERH0[, ww], t='l', lty=4, col=cols)
legend("topright", col=cols, as.character(rhos[ww]), title=expression(rho), lty=1)
legend("topleft", col=1, lty=c(1:2, 4), c("Simes", "kFWER", "kFWER-H0"), horiz=TRUE)
dev.off()

if (FALSE) {
    matplot(-log10(resSimes[, ww]), t='l', lty=1)
    matlines(-log10(resKFWER[, ww]), t='l', lty=2)
    matlines(-log10(resKFWERH0[, ww]), t='l', lty=4)
}

## }} ## end 'for:s'

## TODO: plots on log-scale
