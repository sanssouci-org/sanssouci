library("tidyverse")
library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
#plan(remote, workers = rep("bernoulli", 100))
#plan(eager)

plan(multiprocess, workers = 100)

res <- listenv()
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
    
    if (!file.exists(mpathname)) {
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
                testStepDown(m, dep, B, pi0, SNR, "constant", alphas, kMaxs=kMaxs, flavor="equi-perm")
            }
        }
        ## summarize into JFWER and Power estimates
        mres <- reshape2::melt(as.list(res))  ## how to do this with purrr::map or friends?
        names(mres) <- c("risk", "flavor", "value", "family", "kMax", "alpha", "sid")
        saveRDS(mres, file=mpathname)
        ## NB: 'mres' is not tidy!
    } else {
        mres <- readRDS(mpathname)
        ## summarize into JFWER and Power estimates
        cres <- reshape2::dcast(mres, kMax+alpha+family+flavor~risk, mean, value.var="value", na.rm=TRUE)
        dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, cres)
        head(dat)
        saveRDS(dat, file=pathname)
        
        if (FALSE) {
            JR <- mres %>% 
                spread(risk, value) %>%
                group_by(kMax, alpha, family, flavor) %>%
                summarise(JR = mean(JR, na.rm=TRUE))
            dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, as.data.frame(JR))
            head(dat)
        }
        
        ## stratified Monte-Carlo
        sJR <- mres %>% 
            mutate(alpha=as.numeric(alpha)) %>% 
            filter(risk=="JR") %>%
            spread(flavor, value) %>%
            group_by(kMax, alpha, family) %>%
            mutate(weight=alpha*(Oracle==0) + (1-alpha)*(Oracle==1)) %>%
            mutate(id=Oracle) %>%
            gather(flavor, value, -risk, -family, -kMax, -alpha, -sid, -weight, -id) %>%
            spread(risk, value) %>%
            group_by(kMax, alpha, family, flavor) %>%    
            summarise(
                aa = unique(alpha),  ## of length 1 as we are grouping by alpha
                sJR = (1-aa)*sum((JR==1)*(id==0))/sum(id==0) + 
                    aa*sum((JR==1)*(id==1))/sum(id==1),
                JR = mean(JR, na.rm=TRUE)) %>%
            select(-aa)
        ## frankly, this is getting fairly obscure! probably not the best tidy-way
        
        dat <- cbind(pi0=pi0, dep=dep, SNR=SNR, as.data.frame(sJR))
        head(dat)
        saveRDS(dat, file=spathname)
    }
}

## tidy results and save to single file
fls <- list.files(path, full.names=TRUE)
fls <- fls[-grep("(molten|SMC)", fls)]
id <- gsub(".rds$", "", basename(fls))
names(fls) <- id
dat <- plyr::ldply(fls, readRDS, .id="id")
names(dat)[match("dep", names(dat))] <- "rho"

head(dat)
filename <- sprintf("%s,%s.rds", sname, pname)
pathname <- file.path(resPath, filename)
saveRDS(dat, file=pathname)
print(pathname)

## tidy results and save to single file (stratified MC)
fls <- list.files(path, full.names=TRUE)
fls <- fls[grep("SMC\\.rds$", fls)]
id <- gsub(".rds$", "", basename(fls))
names(fls) <- id
dat <- plyr::ldply(fls, readRDS, .id="id")
names(dat)[match("dep", names(dat))] <- "rho"

head(dat)
filename <- sprintf("%s,%s,SMC.rds", sname, pname)
pathname <- file.path(resPath, filename)
saveRDS(dat, file=pathname)
print(pathname)
