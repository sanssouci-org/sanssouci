rootPath <- "resData/stepDownEqui+power,Mein2006,sansSouci_0.4.1,mine"
##rootPath <- "resData/stepDownEqui+power,equi,sansSouci_0.4.1,rei1"
##rootPath <- "resData/stepDownEqui+power,sansSouci_0.4.1,rei1"

figName <- gsub(".*,(.*),sansSouci.*,.*", "\\1", rootPath)
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

dat <- plyr::ldply(fls, readRDS, .id="id")
dat$se <- sqrt(dat$value * (1-dat$value))/sqrt(nbSimu)

risks <- c("JFWER", "Power")
names(risks) <- risks
datList <- lapply(risks, function(rr) subset(dat, risk==rr))

library("ggplot2")

methods <- c("Simes", "kFWER")
for (rr in risks) {
    datRR <- datList[[rr]]
    str(datRR)

    for (mm in methods) {
        filename <- sprintf("%s,%s,%s,%s.pdf", figName, pname, rr, mm)
        pathname <- file.path(figPath, filename)

        datMM <- datRR[grep(mm, datRR$method), ]

        pdf(pathname)
        p <- ggplot(datMM, aes(x=SNR, y=value, group=method, color=method))
        p <- p + geom_line() + geom_point()
        p <- p + geom_errorbar(aes(ymax=value+se, ymin=value-se), width=0.1)
        p <- p + facet_grid(pi0 ~ rho)
        ##p <- p + facet_grid(pi0 ~ rho, labeller=labeller(rho=label_bquote(rho*"="*.(x)), pi0=label_bquote(pi[0]*"="*.(x))))
        ##    p <- p + scale_x_continuous(breaks=unique(d$SNR)) + xlab(expression(mu))
        p <- p + ylab(rr)
        if (rr=="JFWER") {
            p <- p + scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.25), limits=c(0, 0.3))
            p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
            p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
        } else {
            p <- p + ylim(c(0, 1))
        }

        print(p)
        dev.off()
    }
}

