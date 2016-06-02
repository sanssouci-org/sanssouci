##library("sansSouci")
library("actuar")
library("parallel")
ncores <- 36

if (Rcpp) {
    sname0 <- "stepDownEqui+power,Rcpp"
} else {
    sname0 <- "stepDownEqui+power"
}

ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
htag <- system("hostname", inter=TRUE)
sname <- sprintf("%s,%s,%s", sname0, ptag, htag)

for (SNR in SNRs) {
    for (rho in rhos) {
        for (pi0 in pi0s) {

rootPath <- "resData"

if (FALSE) {
    tags <- sprintf("m=%s,pi0=%s,rho=%s,SNR=%s,B=%s,alpha=%s,nbSimu=%s",
                    m, pi0, rho, SNR, B, alpha, nbSimu)
    filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
    path <- file.path(rootPath, sname)
    path <- R.utils::Arguments$getWritablePath(path)
} else {
    ## one folder per simulation setting, one file per simulation run
    tags <- sprintf("m=%s,pi0=%s,rho=%s,SNR=%s,B=%s,alpha=%s",
                    m, pi0, rho, SNR, B, alpha)
    simdir <- sprintf("%s,nbSimu=%s", gsub("\\.", "_", tags), nbSimu)
#    path <- file.path(rootPath, sname, simdir)
#    path <- R.utils::Arguments$getWritablePath(path)
}

filename <- sprintf("%s,%s.rds", sname, gsub("\\.", "_", tags))
path <- "resData"
path <- file.path(path, sname)
path <- R.utils::Arguments$getWritablePath(path)

print(tags)

res <- mclapply(1:nbSimu, FUN=function(bb) {
    if (bb %% 100 == 0) {
        print(bb)
    }
    
    rr <- testStepDown(m, rho, B, pi0, SNR, alpha)
    return(rr)
}, mc.cores=ncores)

pathname <- file.path(path, filename)
saveRDS(res, file=pathname)

        }
    }
}



