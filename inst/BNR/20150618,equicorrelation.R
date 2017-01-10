library("sansSouci")

## model parameters
m <- 5e2+1
B <- 5e3

## JFWER setup
maxSteps <- 100
methods <- c("Simes", "kFWER")
kMax <- m
alpha <- 0.2

flavor <- "equi-correlated"

rhos <- seq(from=0, to=0.9, by=0.1)

nbSimu <- 10

library("parallel")
if (FALSE) {
    if (detectCores()<30) stop("This is intended to be run on a 32 CPU machine")
    resP <- mclapply(1:length(rhos), FUN=function(rr) {
                         rho <- rhos[rr]
                         print(rho)

                         resRR <- mclapply(1:nbSimu, FUN=function(bb) {
                                               mat <- simulateGaussianNullsFromFactorModel(m, n=B, flavor=flavor, rho=rho)
                                               gammas <- numeric(0L)
                                               for (mm in methods) {
                                                   resJ <- getJointFWERThresholds(mat, tau=mm, alpha, maxSteps=maxSteps, kMax=kMax)
                                                   gammas <- c(gammas, resJ$lambda)
                                               }
                                               gammas
                                           }, mc.cores=8)
                         resRR <- do.call(cbind, resRR)
                     }, mc.cores=3)
    library("abind")
    res <- abind(resP, along=0)
    dimnames(res) <- list(rho=as.character(rhos),
                          method=methods,
                          B=as.character(1:nbSimu))

    saveRDS(res, file="resData/gamma.rds")
}


##
res <- readRDS("resData/gamma.rds")

pdf("../../../../fig/BNR/gamma,Simes.pdf")
par(mar=c(4, 5, 0, 0)+0.1, cex.lab=1.5)
resS <- t(res[, "Simes",])
boxplot(resS, xlab=expression(rho), ylab=expression(lambda(alpha)/alpha))
dev.off()

pdf("../../../../fig/BNR/gamma,kFWER.pdf")
par(mar=c(4, 5, 0, 0)+0.1, cex.lab=1.5)
resK <- t(res[, "kFWER",])
boxplot(resK, xlab=expression(rho), ylab=expression(lambda(alpha)/alpha), ylim=c(0, max(resK)))
abline(h=1/m, lty=2)
dev.off()

## unbalancedness
seqk <- 1:m
alpha <- 0.1

pSimes <- function(k, m, alpha) pbeta(alpha*k/m, k, m-k+1)
pkSimes <- function(alpha, m) sapply(seqk, pSimes, m, alpha)
probakSimes <- pkSimes(alpha, m)
                                        #plot(seqk[1:10],log(probakSimes)[1:10],type="l",col="red",lwd=2,ylab="",xlab="")
ks <- c(1:5,10,100)
cors <- c(1:6)
signif(probakSimes[ks], 2)

lambdas <- alpha*rowMeans(res[, "Simes",])  ## average value of $\lambda(\alpha)$
las <- sapply(lambdas, pkSimes, m)
signif(las[ks, cors], 2)

## kfwer
resK <- mclapply(1:length(rhos), FUN=function(rr) {
    rho <- rhos[rr]
    print(rho)

    mat <- simulateGaussianNullsFromFactorModel(m, n=B, flavor=flavor, rho=rho)

    resList <- NULL
    for (mm in methods) {
        resJ <- getJointFWERThresholds(mat, tau=mm, alpha, maxSteps=maxSteps, kMax=kMax)
        resList[[mm]] <- resJ$probk
    }
    resList
}, mc.cores=3)

library("abind")
resSimes <- sapply(resK, "[[", "Simes")
resKFWER <- sapply(resK, "[[", "kFWER")
colnames(resSimes) <- rhos
colnames(resKFWER) <- rhos

ww <- match(c(0, 0.1, 0.2, 0.5), rhos)
cols <- seq(along=ww)

pdf("unbalanced.pdf")
par(lwd=1.5)
matplot(resSimes[, ww], t='l', lty=1, col=cols, ylab=expression(P(p[(k)] < tau[k])), xlab="k")
matlines(resKFWER[, ww], t='l', lty=2, col=cols)
legend("topright", col=cols, as.character(rhos[ww]), title=expression(rho), lty=1)
legend("topleft", col=1, lty=1:2, c("Simes", "kFWER"))
dev.off()

pdf("unbalanced,zoom.pdf")
matplot(resSimes[, ww], t='l', lty=1, col=cols, ylab=expression(P(p[(k)] < tau[k])), xlim=c(1,20), xlab="k")
matlines(resKFWER[, ww], t='l', lty=2, col=cols)
legend("topright", col=cols, as.character(rhos[ww]), title=expression(rho), lty=1)
legend("topleft", col=1, lty=1:2, c("Simes", "kFWER"))
dev.off()

matplot(-log10(resSimes[, ww]), t='l', lty=1)
matlines(-log10(resKFWER[, ww]), t='l', lty=2)

