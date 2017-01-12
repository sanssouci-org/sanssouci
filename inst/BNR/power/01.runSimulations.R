library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(remote, workers = rep("bernoulli", 100))

for (ii in 1:nrow(configs)) {
    pi0 <- configs[ii, "pi0"]
    dep <- configs[ii, "dep"]
    kMaxs <- c(max(round(m*(1-pi0)), 1), m/2, m)
    
    tags <- sprintf("pi0=%s,dep=%s", pi0, dep)
    print(tags)
    filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
    pathname <- file.path(path, filename)
    if (file.exists(pathname)) {
        next;
    }    
    
    ## using futures =)
    res <- listenv()
    options("future.wait.times"=1e4)
    for (ss in 1:nbSimu) {
        if (ss %% 100==0) { print(ss);}
        res[[ss]] %<-% {
            library("sansSouci")
            pathname <- system.file("BNR/testStepDown.R", package="sansSouci")
            stopifnot(file.exists(pathname))
            source(pathname); rm(pathname)
            testStepDown(m, dep, B, pi0, SNR, typeOfSNR, alphas, kMaxs=kMaxs, flavor=flavor)
        }
    }
    print(head(res))
    ## summarize into JFWER and Power estimates
    mres <- reshape2::melt(as.list(res))
    names(mres) <- c("risk", "flavor", "value", "family", "kMax", "alpha", "sid")
    
    ## summarize into JFWER and Power estimates
    cres <- reshape2::dcast(mres, kMax+alpha+family+flavor~risk, mean, value.var="value")
    dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, cres)
    head(dat)
    
    saveRDS(dat, file=pathname)
}
