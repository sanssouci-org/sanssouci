##library("sansSouci")
#library("actuar")
library("parallel")
ncores <- 38

ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
htag <- system("hostname", inter=TRUE)
sname <- sprintf("%s,%s,%s", sname0, ptag, htag)

kMax <- m
SNR <- 2
pname <- sprintf("m=%s,B=%s,SNR=%s,nbSimu=%s", m, B, SNR, nbSimu)

resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

alphas <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
for (alpha in alphas) {
    for (dep in deps) {
        for (pi0 in pi0s) {
            tags <- sprintf("pi0=%s,kMax=%s,dep=%s,SNR=%s", pi0, kMax, dep, SNR)
            filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
            print(tags)

            res <- mclapply(1:nbSimu, FUN=function(bb) {
                if (bb %% 100 == 0) {
                    print(bb)
                }

                rr <- testStepDown(m, dep, B, pi0, SNR, typeOfSNR, alpha, kMax=kMax, flavor=flavor)
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
                eV0 <- colMeans(dat[which(nms=="v0"), ])
                eS1 <- colMeans(dat[which(nms=="s1"), ])

                ## TODO: add sds, ie
                if (FALSE) {
                    d <- dat[which(nms=="rej0"), ]>0
                    mode(d) <- "numeric"
                    matrixStats::colSds(d)
                }

                mat <- rbind(JFWER=eJR, Power=ePow, V0=eV0, S1=eS1)
                names(dimnames(mat)) <- c("risk", "method")
                dat <- reshape2::melt(mat)
                dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, dat)
            }
            pathname <- file.path(path, filename)
            saveRDS(dat, file=pathname)
        }
    }
}
