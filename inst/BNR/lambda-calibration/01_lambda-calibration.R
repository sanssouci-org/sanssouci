oneSim <- function(m, rho, B, alpha, pi0=1, SNR=0, kMaxs=kMaxs) {
    sim <- gaussianTestStatistics(m, B, pi0 = pi0, SNR = SNR, dep = "equi", param = rho)
    configs <- expand.grid(kMax=kMaxs, refFam=c("kFWER", "Simes"), stringsAsFactors=FALSE)
    lambdas <- numeric(0L)
    for (rr in 1:nrow(configs)) {
        fam <- configs[rr, "refFam"]
        kMax <- configs[rr, "kMax"]
        resJ <- jointFWERControl(sim$X0, refFamily=fam, alpha, kMax=kMax, verbose=FALSE)
        lambdas <- c(lambdas, resJ$lambda)
    }
    df <- cbind(configs, lambda=lambdas)
    return(df)
}

resr <- listenv()
for (rr in seq(along = rhos)) {
    rho <- rhos[rr]
    res <- listenv()
    for (ss in 1:nbSimu) {
        if (ss %% 20==0) { print(ss);}
        res[[ss]] %<-% {
            library("sansSouci")
            oneSim(m, rho, B, alpha, kMaxs=kMaxs)
        }
    }
    names(res) <- 1:nbSimu
    dat <- plyr::ldply(as.list(res), data.frame, .id="sid")
    resr[[rr]] <- dat
}
names(resr) <- rhos
lambdas <- plyr::ldply(resr, data.frame, .id="rho")

#saveRDS(lambdas, file="resData/lambdas.rds")
