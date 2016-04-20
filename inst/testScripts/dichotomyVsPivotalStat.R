m <- 10230
B <- 1e4

flavor <- c("independent", "equi-correlated", "3-factor model")[2]
rho <- 0.2

mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
str(mat)

alpha <- 0.2

maxSteps <- 100
kMax <- m

refFamily <- "Simes"
flavor <- "dichotomy"
system.time(res2 <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                           maxSteps=maxSteps, kMax=kMax,
                                           flavor=flavor))
flavor <- "pivotalStat"
system.time(res <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                          maxSteps=maxSteps, kMax=kMax,
                                          flavor=flavor))

refFamily <- "kFWER"
flavor <- "dichotomy"
system.time(resk2 <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                           maxSteps=maxSteps, kMax=kMax,
                                           flavor=flavor))
flavor <- "pivotalStat"
system.time(resk <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                           maxSteps=maxSteps, kMax=kMax))


refFamily <- "kFWER"
flavor <- "dichotomy"
system.time(resk2c <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                           maxSteps=maxSteps, kMax=kMax,
                                           flavor=flavor, Rcpp=TRUE))



##
foo <- function(m, B, alpha) {
    mat <- simulateGaussianNullsFromFactorModel(m, B, flavor="independent")
    str(mat)
    system.time(ps <- getJointFWERThresholds(mat, refFamily="kFWER", alpha, flavor="pivotalStat"))
    system.time(di <- getJointFWERThresholds(mat, refFamily="kFWER", alpha))
    plot(sort(ps$pivotalStat)[round(alpha*B)-(-100:100)])
    abline(h=c(ps$lambda, di$lambda), col=c(1,2), lty=c(1,4))
    list(ps, di)
}

#devtools::install_github("rstudio/profvis")
library("profvis")

profvis({
    resk2c <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                     maxSteps=maxSteps, kMax=kMax,
                                     flavor=flavor, Rcpp=TRUE)
})

p <- profvis({
    resk2 <- getJointFWERThresholds(mat, refFamily=refFamily, alpha,
                                     maxSteps=maxSteps, kMax=kMax,
                                     flavor=flavor, Rcpp=FALSE)
})
