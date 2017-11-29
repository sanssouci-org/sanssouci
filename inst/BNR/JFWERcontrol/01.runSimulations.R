library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
#plan(remote, workers = rep("bernoulli", 100))
#plan(eager)

plan(multiprocess, workers = 100)
#res <- listenv()
for (ii in 1:nrow(configs)) {
    pi0 <- configs[ii, "pi0"]
    dep <- configs[ii, "dep"]
    SNR <- configs[ii, "SNR"]
    
    kMaxs <- c(max(2, round(2*m*(1-pi0))), 2, 10, m)
    kMaxs <- unique(kMaxs)
    
    tags <- sprintf("pi0=%s,dep=%s,SNR=%s", pi0, dep, SNR)
    print(tags)
    filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
    pathname <- file.path(path, filename)
    sfilename <- sprintf("%s,SMC.rds", gsub("\\.", "_", tags)) ## stratified Monte-Carlo
    spathname <- file.path(path, sfilename)
    mfilename <- sprintf("%s,molten.rds", gsub("\\.", "_", tags))
    mpathname <- file.path(path, mfilename)
    
    if (file.exists(mpathname)) {
        next;
    } else {
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
                testStepDown(m, dep, B, pi0, SNR, "constant", alphas, kMaxs=kMaxs, flavor=simFlavor)
            }
        }
        ## summarize into JFWER and Power estimates
        mres <- reshape2::melt(as.list(res))  ## how to do this with purrr::map or friends?
        names(mres) <- c("risk", "flavor", "value", "family", "kMax", "alpha", "sid")
        saveRDS(mres, file=mpathname)
        ## NB: 'mres' is not tidy!
    } 
}

