library("sansSouci")

typeOfSNR <- "constant"
flavor <- "equi"
sname0 <- sprintf("%s,%s", flavor, typeOfSNR)
ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
sname <- sprintf("%s,%s", sname0, ptag)

SNRs <- rev(c(0, 1, 2, 3, 4, 5))
deps <- c(0, 0.2, 0.4) ## correlation coefficient
pi0s <- c(0.8, 0.9, 0.99, 0.999)
pi0s <- c(0.8, 0.9, 0.99)
alphas <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)

m <- 1e3
B <- 1e3
nbSimu <- 1e3


library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(cluster, workers = rep("bernoulli", 101))

kMax <- m
SNR <- 2
pname <- sprintf("m=%s,B=%s,SNR=%s,nbSimu=%s", m, B, SNR, nbSimu)

resPath <- "resData"
path <- file.path(resPath, sname, pname)
path <- R.utils::Arguments$getWritablePath(path)

##configs <- expand.grid(alpha=alphas, pi0=pi0s, dep=deps)
configs <- expand.grid(pi0=pi0s, dep=deps)
configs

for (ii in 1:nrow(configs)) {
    pi0 <- configs[ii, "pi0"]
    dep <- configs[ii, "dep"]
    kMaxs <- c(max(round(m*(1-pi0)), 1), m/2, m)
    
    tags <- sprintf("pi0=%s,dep=%s", pi0, dep)
    filename <- sprintf("%s.rds", gsub("\\.", "_", tags))
    print(tags)

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
    
    pathname <- file.path(path, filename)
    saveRDS(dat, file=pathname)
}


## plot
fls <- list.files(path, full.names=TRUE)
id <- gsub(".rds$", "", basename(fls))
names(fls) <- id
dat <- plyr::ldply(fls, readRDS, .id="id")

risks <- c("JFWER", "Power", "v0", "s1")
risks <- c("JR", "Power", "v0", "s1")

figName <- sname0
figPath <- file.path("fig", sprintf("BNR,%s", ptag))
figPath <- R.utils::Arguments$getWritablePath(figPath)

library("ggplot2")

methods <- c("Simes", "kFWER")
for (rr in risks) {
    for (mm in methods) {
        filename <- sprintf("%s,%s,%s,%s.pdf", figName, pname, rr, mm)
        pathname <- file.path(figPath, filename)

        datMM <- dat[grep(mm, dat$family), ]

        pdf(pathname)
        p <- ggplot(datMM, aes(x=SNR, y=eval(rr), group=family, color=family))
        p <- p + geom_line() + geom_point()
        ##p <- p + geom_errorbar(aes(ymax=value+se, ymin=value-se), width=0.1)
        p <- p + facet_grid(pi0 ~ rho)
        ##p <- p + facet_grid(pi0 ~ rho, labeller=labeller(rho=label_bquote(rho*"="*.(x)), pi0=label_bquote(pi[0]*"="*.(x))))
        ##    p <- p + scale_x_continuous(breaks=unique(d$SNR)) + xlab(expression(mu))
        p <- p + ylab(rr)
        if (rr=="JR") {
            ##            p <- p + scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.25), limits=c(0, 0.35))
            p <- p + scale_y_continuous(breaks=c(0, alpha), limits=c(0, alpha*1.3))
            p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
            p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
        } else {
            p <- p + ylim(c(0, 1))
        }

        print(p)
        dev.off()
    }
}

