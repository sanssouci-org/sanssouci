library("sansSouci")

sname0 <- "JFWERcontrol"
ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
sname <- sprintf("%s,%s", sname0, ptag)

deps <- c(0, 0.2, 0.4) ## correlation coefficient
pi0s <- c(0.8, 0.9, 0.99, 0.999)
SNRs <- c(0, 1, 2, 3, 4, 5)[c(1, 3, 6)]

alphas <- c(0.25)

m <- 1e3
B <- 1e3
nbSimu <- 1e3

pname <- sprintf("m=%s,B=%s,nbSimu=%s", m, B, nbSimu)
resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

configs <- expand.grid(pi0=pi0s, dep=deps, SNR=SNRs)
configs

