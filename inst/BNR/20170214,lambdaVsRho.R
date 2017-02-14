## For Figures 2 and 3 in the BNR paper
library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")

library("sansSouci")

## setup
m <- 1e3
kMaxs <- c(m, 10)
B <- 1e4
methods <- c("Simes", "kFWER")
alpha <- 0.2

nbSimu <- 100
rhos <- seq(from=0, to=1, by=0.1)

oneSim <- function(m, rho, B, alpha, pi0=1, SNR=0, kMaxs=kMaxs) {
    sim <- sansSouci::simulateEqui(m, rho, B, pi0=pi0, SNR=SNR)
    configs <- expand.grid(kMax=kMaxs, refFam=c("kFWER", "Simes"), stringsAsFactors=FALSE)
    lambdas <- numeric(0L)
    for (rr in 1:nrow(configs)) {
        fam <- configs[rr, "refFam"]
        kMax <- configs[rr, "kMax"]
        resJ <- jointFWERControl(sim$X0, refFamily=fam, alpha, kMax=kMax, verbose=FALSE)
        lambdas <- c(lambdas, resJ$lambda)
    }
    df <- cbind(configs, lambda=lambdas)
    return(df)
}

#plan(remote, workers = rep("bernoulli", 100))
plan(multiprocess, workers = 100)


res <- listenv()
for (rho in rhos) {
    options("future.wait.times"=1e4)
    for (ss in 1:nbSimu) {
        if (ss %% 20==0) { print(ss);}
        res[[ss]] %<-% {
            library("sansSouci")
            oneSim(m, rho, B, alpha, kMaxs=kMaxs)
        }
    }
}
rdat <- Reduce(rbind, res)
rdat <- dplyr::bind_rows(as.list(res), .id=1:10)
stop("Done!")



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


