library("sansSouci")

## model parameters
m <- 1e3
B <- 1e4

## JFWER setup
maxSteps <- 100
methods <- c("Simes", "kFWER")
kMax <- m
alpha <- 0.2

flavor <- "equi-correlated"

rhos <- seq(from=0, to=0.9, by=0.1)

nbSimu <- 10

library("parallel")
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
    res[rr,,] <- resRR
}, mc.cores=3)

library("abind")
res <- abind(resP, along=0)
dimnames(res) <- list(rho=as.character(rhos),
                      method=methods,
                      B=as.character(1:nbSimu))

saveRDS(res, file="resData/gamma.rds")
res <- readRDS("resData/gamma.rds")
boxplot(t(res[, 1,]), ylab=expression(gamma), xlab=expression(rho), main="gamma adjustment -- Simes thresholds")

