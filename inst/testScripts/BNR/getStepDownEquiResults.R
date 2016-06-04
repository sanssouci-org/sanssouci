rootPath <- "resData/stepDownEqui+power,Rcpp,sansSouci_0.3.1,rei1"
pname <- sprintf("m=%s,B=%s,alpha=%s,nbSimu=%s", m, B, alpha, nbSimu)
path <- file.path(rootPath, pname)
path <- R.utils::Arguments$getReadablePath(path)

library("R.utils")
figPath <- "fig/BNR"
figPath <- Arguments$getWritablePath(figPath)

## pattern <- sprintf("stepDownEqui+power,Rcpp,m=%s,pi0=%s,rho=%s,SNR=%s,B=%s,alpha=%s,nbSimu=%s",
##                    m, pi0, rho, SNR, B, alpha, nbSimu)


fls <- list.files(path, full.names=TRUE)
id <- gsub(".rds$", "", basename(fls))
names(fls) <- id
str(fls)
patt <- "pi0=(.*),rho=(.*),SNR=([0-9]+)"

if (FALSE) {  ## for a previous version of the saved results
    resList <- lapply(fls, FUN=function(ff) {
        res <- readRDS(ff)
        dat <- Reduce(rbind, res)
        nms <- rownames(dat)
        eJR <- colMeans(dat[which(nms=="rej0"), ]>0)
        ePow <- colMeans(dat[which(nms=="rej1"), ]>0)
        mat <- rbind(JFWER=eJR, Power=ePow)
        names(dimnames(mat)) <- c("risk", "method")
        dat <- reshape2::melt(mat)
        
        id <- gsub(".rds$", "", basename(ff))
        pi0 <- gsub(patt, "\\1", id)
        rho <- gsub(patt, "\\2", id)
        SNR <- gsub(patt, "\\3", id)
        res <- cbind(pi0=pi0, rho=rho, SNR=SNR, dat)
        saveRDS(res, file=ff)
    })
}

dat <- plyr::ldply(fls, readRDS, .id="id")
dat$se <- sqrt(dat$value * (1-dat$value))/sqrt(nbSimu)

risks <- c("JFWER", "Power")
names(risks) <- risks
datList <- lapply(risks, function(rr) subset(dat, risk==rr))

figName <- "stepDownEqui"

for (rr in risks) {
    filename <- sprintf("%s,%s,%s.pdf", figName, pname, rr)
    pathname <- file.path(figPath, filename)
    datRR <- datList[[rr]]
    str(datRR)
    
    pdf(pathname)
    p <- ggplot(datRR, aes(x=SNR, y=value, group=method, color=method))
    p <- p + geom_line() + geom_point() 
    p <- p + geom_errorbar(aes(ymax=value+se, ymin=value-se), width=0.1)
    p <- p + facet_grid(pi0 ~ rho)
    ##p <- p + facet_grid(pi0 ~ rho, labeller=labeller(rho=label_bquote(rho*"="*.(x)), pi0=label_bquote(pi[0]*"="*.(x))))
    ##    p <- p + scale_x_continuous(breaks=unique(d$SNR)) + xlab(expression(mu))
    p <- p + scale_y_continuous(breaks=c(0.1, 0.2, 0.25))
    p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
    p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
    print(p)
    dev.off()
}

