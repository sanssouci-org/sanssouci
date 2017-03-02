library("sansSouci")

sname0 <- "detectionPower"
ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
sname <- sprintf("%s,%s", sname0, ptag)

deps <- c(0, 0.2, 0.4) ## correlation coefficient
alphas <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
# betar <- rbind(
#     cbind(beta=0.5, r=c(0.1, 0.2)),
#     cbind(beta=2/3, r=c(0.1, 0.2, 0.3, 0.4, 0.5)),
#     cbind(beta=0.8, r=c(0.1, 0.2, 0.3, 0.4, 0.5)),
#     cbind(beta=1.0, r=c(0.1, 0.5, 1, 2)))
# pi1 <- m^(-betar[, "beta"])
# signif(unique(pi1), 1)

betas <- c(0.5, 0.6, 2/3, 0.8, 1)
rs <- c(0.05, 0.1, 0.2, 0.5, 1)

m <- 1e3
B <- 1e3
nbSimu <- 1e3

pname <- sprintf("m=%s,B=%s,nbSimu=%s", m, B, nbSimu)
resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

pi1 <- m^(-betas)
signif(unique(pi1), 1)

configs <- expand.grid(beta=betas, r=rs, dep=deps)
configs
