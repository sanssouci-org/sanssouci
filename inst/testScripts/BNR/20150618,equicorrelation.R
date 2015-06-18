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

rhos <- c(0, 0.1, 0.2, 0.3, 0.5, 0.8)
parName <- "rho"
params <- rhos

lambdas <- NULL
probs <- NULL

nbSimu <- 6
res <- array(NA_real_, dim=c(length(rhos), length(methods), nbSimu))
dimnames(res) <- list(rho=as.character(rhos),
                      method=methods,
                      B=as.character(1:nbSimu))

library("parallel")
for (rr in seq(along=rhos)) {
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
                      }, mc.cores=3)
    resRR <- do.call(cbind, resRR)
    res[rr,,] <- resRR
}



