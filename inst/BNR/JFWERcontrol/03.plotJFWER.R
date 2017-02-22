filename <- sprintf("%s,%s.rds", sname, pname)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
head(dat)

gat <- tidyr::gather(dat, "criterion", "value", JR, detPow, detPow1, v0, estPow, estPow1, powBH5, powBH50, pow0)

## some reshaping
datC <- subset(gat, criterion=="JR")
datC$family <- factor(datC$family, levels=c("kFWER", "Simes"), labels=c("Balanced", "Linear"))
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)

kc <- as.character(datC$kMax)
kc[which(datC$kMax==m)] <- "m"
kc[which(datC$kMax==m/2)] <- "m/2"
kc[which(datC$kMax==2*(1-datC$pi0)*m)] <- "2m1"
datC$kMaxC <- kc

datC$ff <- datC$flavor

figName <- sname0
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names

## plot various statistics vs nominal JFWER
library("ggplot2")

kMaxs <- unique(datC$kMaxC)
families <- unique(datC$family)
confs <- expand.grid(alpha=alphas, kMax=kMaxs, family=families, stringsAsFactors=FALSE)

for (ii in 1:nrow(confs)) {
    kk <- confs[ii, "kMax"]
    aa <- confs[ii, "alpha"]
    fam <- confs[ii, "family"]
    ftag <- sprintf("family=%s,alpha=%s,kMax=%s", fam, aa, kk)
    
    filename <- sprintf("%s,BalancedVsLinear,%s,%s.pdf", figName, pname2, ftag)
    pathname <- file.path(figPath, filename)
    datI <- subset(datC, family==fam & alpha==aa & kMaxC==kk & flavor != "Oracle2")
    pal <- RColorBrewer::brewer.pal(8,"BrBG")
    
    if (fam=="Balanced") {
        datI <- subset(datI, flavor != "unadjusted")
        pal <- pal[1:3]
    } else if (fam=="Linear") {
        pal <- pal[5:8]
        datI$ff <- plyr::revalue(datI$ff, c("unadjusted"="Simes"))
    }
    pdf(pathname)
    ##
    p <- ggplot(datI, aes_string(x="SNR", y="value", group="ff", color="ff"))
    p <- p + geom_line()
    p <- p + facet_grid(pi0 ~ rho,
                        #    scales="free_y",
                        labeller=label_bquote(
                            rows= pi[0]==.(pi0),
                            cols= rho==.(rho)))
    p <- p + scale_y_continuous(breaks=c(0, aa), limits=c(0, aa*1.3))
    p <- p + labs(color=sprintf("%s family", fam))
    p <- p + geom_point()
    p <- p + labs(y="JR", x=expression(bar(mu)))
    p <- p + scale_colour_manual(values = rev(pal))
    p <- p + theme(axis.title=element_text(size=16),
                   strip.text = element_text(size=12),
                   legend.text = element_text(size=12),
                   legend.title = element_text(size=12))
    print(p)
    dev.off()
}

