library("sansSouci")
## source("package/R/getJointFWERThresholds.R")
## Rcpp::sourceCpp("package/src/coverage.cpp")
## Rcpp::sourceCpp("package/src/partialColSort.cpp")

m <- 1023
B <- 1e4

flavor <- c("independent", "equi-correlated", "3-factor model")[2]
rho <- 0.2
mat <- simulateGaussianNullsFromFactorModel(m, B, flavor=flavor, rho=rho)
alpha <- 0.2

Rprof("jfwer.Rout")
res <- getJointFWERThresholds(mat, refFamily="kFWER", alpha)
Rprof(NULL)
str(res)

sapply(summaryRprof("jfwer.Rout"), head)


Rprof("jfwer,C.Rout")
resC <- getJointFWERThresholds(mat, refFamily="kFWER", alpha, Rcpp=TRUE)
Rprof(NULL)
str(resC)

sapply(summaryRprof("jfwer,C.Rout"), head)

## same results?
if (FALSE) {
    for (kk in names(res)) {
        print(kk);
        print(identical(res[[kk]], resC[[kk]]))
    }

    kk <- "sLambda"
    alphas <- 1:100/100
    diffs <- sapply(alphas, FUN=function(alpha) {
                        max(abs(resC$sLambda(alpha)-res$sLambda(alpha)))
                    })
    max(abs(diffs))
    ## => identical results
}
