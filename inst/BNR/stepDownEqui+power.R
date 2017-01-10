##library("sansSouci")
#library("actuar")
library("parallel")
library("reshape2")
ncores <- 3

ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
htag <- system("hostname", inter=TRUE)
sname <- sprintf("%s,%s,%s", sname0, ptag, htag)
pname <- sprintf("m=%s,B=%s,alpha=%s,nbSimu=%s", m, B, alpha, nbSimu)

stag <- sprintf("%s,%s,m=%s,B=%s,nbSimu=%s", flavor, typeOfSNR, m, B, nbSimu)
resPath <- sprintf("resData,%s", htag)
path <- file.path(resPath, stag)
path <- R.utils::Arguments$getWritablePath(path)

for (SNR in SNRs) {
    for (dep in deps) {
        for (pi0 in pi0s) {
            kMaxs <- c(round(m*(1-pi0)), m/2, m)
            tags <- sprintf("pi0=%s,dep=%s,SNR=%s", pi0, dep, SNR)
            filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
            print(tags)

            res <- mclapply(1:nbSimu, FUN=function(bb) {
                if (bb %% 10 == 0) {
                    print(bb)
                }

                rr <- testStepDown(m, dep, B, pi0, SNR, typeOfSNR, alphas, kMaxs=kMaxs, flavor=flavor)
            }, mc.cores=ncores)
            mres <- melt(res)
            names(mres) <- c("risk", "flavor", "value", "family", "kMax", "alpha", "sid")

            ## summarize into JFWER and Power estimates
            cres <- dcast(mres, kMax+alpha+family+flavor~risk, mean, value.var="value")
            dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, cres)

            pathname <- file.path(path, filename)
            saveRDS(dat, file=pathname)
        }
    }
}
