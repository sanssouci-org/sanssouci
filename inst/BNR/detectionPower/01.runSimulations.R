library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
#plan(remote, workers = rep("bernoulli", 100))
#plan(eager)
plan(multiprocess, workers = 50)

res <- listenv()
for (ii in 1:nrow(configs)) {
    beta <- configs[ii, "beta"]
    r <- configs[ii, "r"]
    dep <- configs[ii, "dep"]
    pi0 <- 1-m^(-beta)
    
    SNR <- sqrt(2*r*log(m))
    kMaxs <- c(round(2*m*(1-pi0)), 10, 20, m)
    kMaxs <- unique(kMaxs)
    tags <- sprintf("beta=%s,r=%s", beta, r)
    print(tags)
    filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
    pathname <- file.path(path, filename)
    mfilename <- sprintf("%s,molten.rds", gsub("\\.", "_", tags))
    mpathname <- file.path(path, mfilename)
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
            testStepDown(m, dep, B, pi0, SNR, "constant", alphas, kMaxs=kMaxs, flavor=flavor)
        }
    }
    
    ## summarize into JFWER and Power estimates
    mres <- reshape2::melt(as.list(res))
    names(mres) <- c("risk", "flavor", "value", "family", "kMax", "alpha", "sid")
    saveRDS(mres, file=mpathname)
    
    ## summarize into JFWER and Power estimates
    cres <- reshape2::dcast(mres, kMax+alpha+family+flavor~risk, mean, value.var="value", na.rm=TRUE)
    dat <- cbind(beta=beta, r=r, pi0=pi0, dep=dep, SNR=SNR, cres)
    head(dat)
    
    saveRDS(dat, file=pathname)
}
