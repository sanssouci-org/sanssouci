##library("sansSouci")
library("actuar")
library("parallel")
ncores <- 36

ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
htag <- system("hostname", inter=TRUE)
sname <- sprintf("%s,%s,%s", sname0, ptag, htag)
pname <- sprintf("m=%s,B=%s,alpha=%s,nbSimu=%s", m, B, alpha, nbSimu)

resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

for (SNR in SNRs) {
    for (rho in rhos) {
        for (pi0 in pi0s) {
            tags <- sprintf("pi0=%s,rho=%s,SNR=%s", pi0, rho, SNR)
            filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
            print(tags)

            res <- mclapply(1:nbSimu, FUN=function(bb) {
                if (bb %% 100 == 0) {
                    print(bb)
                }

                rr <- testStepDown(m, rho, B, pi0, SNR, alpha)
            }, mc.cores=ncores)

            ## tidy? too many rows! (m*nbSimu>=1e6 per setting...)
            if (FALSE) {
                names(res) <- as.character(1:nbSimu)
                dat <- plyr::ldply(res, cbind, .id="runId")
            } else {  ## summarize into JFWER and Power estimates
                dat <- Reduce(rbind, res)
                nms <- rownames(dat)
                eJR <- colMeans(dat[which(nms=="rej0"), ]>0)
                ePow <- colMeans(dat[which(nms=="rej1"), ]>0)
                mat <- rbind(JFWER=eJR, Power=ePow)
                names(dimnames(mat)) <- c("risk", "method")
                dat <- reshape2::melt(mat)
                dat <- cbind(pi0=pi0, rho=rho, SNR=SNR, dat)
            }
            pathname <- file.path(path, filename)
            saveRDS(dat, file=pathname)
        }
    }
}
