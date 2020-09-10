library("future")
library("sansSouci")
plan(eager)
plan(remote, workers=rep("bernoulli", 1))
res %<-% {
    print("hou")
    example(jointFWERControl)
    
    pathname <- system.file("BNR/00.params.R", package="sansSouci")
    stopifnot(file.exists(pathname))
    source(pathname); rm(pathname)
    alphas <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
    pi0 <- 0.8
    kMaxs <- unique(c(round(m*(1-pi0)), m/2, m))

    testStepDown(1e3, 0, 1e3, 0.8, 1, "constant", alphas, kMaxs=kMaxs)
}
res

     
plan(remote, workers=rep("bernoulli", 20))
library(listenv)
res <- listenv()
for (ss in 1:30) {
    res[[ss]] %<-% {
        pathname <- system.file("BNR/00.params.R", package="sansSouci")
        stopifnot(file.exists(pathname))
        source(pathname); rm(pathname)
        
        testStepDown(1e3, 0, 1e3, 0.8, 1, "constant", 0.1, kMaxs=m)
    }
}
str(res)
