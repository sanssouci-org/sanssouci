library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(remote, workers = rep("bernoulli", 10))

## retrieve results (remote!)
test %<-% getwd()
dat %<-% {
    rpath <- file.path("~/Documents/Packages/sanssouci", path)
    fls <- list.files(rpath, full.names=TRUE)
    id <- gsub(".rds$", "", basename(fls))
    names(fls) <- id
    plyr::ldply(fls, readRDS, .id="id")
}
risks <- c("JFWER", "Power", "v0", "s1")
risks <- c("JR", "Power", "v0", "s1")

## plot
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
        p <- ggplot(datMM, aes(x=SNR, y=JR, group=family, color=family))
        p <- p + geom_line() + geom_point()
        ##p <- p + geom_errorbar(aes(ymax=value+se, ymin=value-se), width=0.1)
        p <- p + facet_grid(pi0 ~ dep)
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

