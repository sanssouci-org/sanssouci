p <- 0.5
m <- 200
nbSimu <- 1e4
B <- 2e3
n <- 100
rho <- 0
pi0 <- 1
alpha <- 0.2

library(howmany.pn)
library(sansSouci)
library(R.utils)
## source("R/getJointFWERThresholds.R")
mpath <- system.file("testScripts/Mein2006/R", package="sansSouci")
sourceDirectory(mpath)

vers <- packageDescription("sansSouci")$Version
flavorH <- "sansSouci"
sort <- TRUE

## paths
study <- ifelse(sort, sprintf("Mein%s,sorted", flavorH), sprintf("Mein%s", flavorH))
if (flavorH=="sansSouci") {
  study <- sprintf("Mein2006,%s,%s", flavorH, vers)
}

path <- file.path("resData", study)
path <- Arguments$getWritablePath(path)

res <- runMeinshausen2006(m, pi0, n, rho, nbSimu, flavorH=flavorH, path=path,
                          overwrite=FALSE, verbose=verbose, sort=sort, alpha=alpha)

