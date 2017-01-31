library("sansSouci")

typeOfSNR <- "constant"
flavor <- "equi"
sname0 <- sprintf("%s,%s", flavor, typeOfSNR)
ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
sname <- sprintf("%s,%s", sname0, ptag)

SNRs <- c(0, 1, 2, 3, 4, 5)
deps <- c(0, 0.2, 0.4) ## correlation coefficient
pi0s <- c(0.8, 0.9, 0.99, 0.999)
alphas <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)

SNRs <- 2:4; deps <- 0

m <- 1e4
B <- 1e3
nbSimu <- 1e3

pname <- sprintf("m=%s,B=%s,nbSimu=%s", m, B, nbSimu)
resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

##configs <- expand.grid(alpha=alphas, pi0=pi0s, dep=deps)
configs <- expand.grid(pi0=pi0s, dep=deps, SNR=SNRs)
configs

