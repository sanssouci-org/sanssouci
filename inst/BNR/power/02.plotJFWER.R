library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(remote, workers = rep("bernoulli", 2)) ## retrieve results (remote!)
plan(eager) ## retrieve results (local!)

test %<-% getwd()
dat %<-% {
    rpath <- file.path("~/Documents/Packages/sanssouci", path)
    fls <- list.files(rpath, full.names=TRUE)
    id <- gsub(".rds$", "", basename(fls))
    names(fls) <- id
    plyr::ldply(fls, readRDS, .id="id")
}
head(dat)

kc <- as.character(dat$kMax)
kc[which(dat$kMax==m)] <- "m"
kc[which(dat$kMax==m/2)] <- "m/2"
kc[which(dat$kMax==(1-dat$pi0)*m)] <- "(1-pi0)*m"
dat$kMaxC <- kc

## plot JFWER
figName <- sname0
figPath <- file.path("fig", sprintf("BNR,%s", ptag))
figPath <- R.utils::Arguments$getWritablePath(figPath)

library("ggplot2")
pal <- c("#DFC27D", "#A6611A", "#018571", "#80CDC1")

alpha <- 0.25
kMax <- m

risks <- "JR"
methods <- c("Simes", "kFWER")
for (rr in risks) {
    for (mm in methods) {
        filename <- sprintf("%s,%s,%s,%s.pdf", figName, pname, rr, mm)
        pathname <- file.path(figPath, filename)
        
        ww <- which(dat$family==mm & dat$alpha==as.character(alpha) & dat$kMax==as.character(kMax))
        datMM <- dat[ww, ]

        pdf(pathname)
        p <- ggplot(datMM, aes(x=SNR, y=JR, group=flavor, color=flavor))
        p <- p + geom_line() + geom_point()
        ##p <- p + geom_errorbar(aes(ymax=value+se, ymin=value-se), width=0.1)
        p <- p + facet_grid(pi0 ~ dep)
        ##p <- p + facet_grid(pi0 ~ rho, labeller=labeller(rho=label_bquote(rho*"="*.(x)), pi0=label_bquote(pi[0]*"="*.(x))))
        ##    p <- p + scale_x_continuous(breaks=unique(d$SNR)) + xlab(expression(mu))
        p <- p + ylab(rr)
        if (rr=="JR") {
            ##            p <- p + scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.25), limits=c(0, 0.35))
            p <- p + scale_y_continuous(breaks=c(0, alpha), limits=c(0, alpha*1.3))
            #p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
            # p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
        } else {
            p <- p + ylim(c(0, 1))
        }
        
        print(p)
        dev.off()
    }
}

