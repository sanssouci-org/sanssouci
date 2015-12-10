m <- 10230
B <- 1e4

flavor <- c("independent", "equi-correlated", "3-factor model")[2]
rho <- 0.2

mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
str(mat)

alpha <- 0.2

maxSteps <- 100
kMax <- m

meth <- "Simes"
flavor <- "dichotomy"
system.time(res2 <- getJointFWERThresholds(mat, tau=meth, alpha,
                                           maxSteps=maxSteps, kMax=kMax,
                                           flavor=flavor))
flavor <- "pivotalStat"
system.time(res <- getJointFWERThresholds(mat, tau=meth, alpha,
                                          maxSteps=maxSteps, kMax=kMax,
                                          flavor=flavor))

meth <- "kFWER"
flavor <- "dichotomy"
system.time(resk2 <- getJointFWERThresholds(mat, tau=meth, alpha,
                                           maxSteps=maxSteps, kMax=kMax,
                                           flavor=flavor))
flavor <- "pivotalStat"
system.time(resk <- getJointFWERThresholds(mat, tau=meth, alpha,
                                          maxSteps=maxSteps, kMax=kMax,



