filename <- sprintf("%s,%s,SMC.rds", sname, pname)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
head(dat)

#gat <- tidyr::gather(dat, "criterion", "value", JR, detPow, detPow1, v0, estPow, estPow1, powBH5, powBH50, pow0)

stratif <- FALSE

## some reshaping
if (stratif) {
    gat <- tidyr::gather(dat, "criterion", "value", "sJR")
    datC <- subset(gat, criterion=="sJR" & pi0 < 0.999)
} else {
    gat <- tidyr::gather(dat, "criterion", "value", "JR")
    datC <- subset(gat, criterion=="JR" & pi0 < 0.999)
}
datC$family <- factor(datC$family, 
                      levels=c("kFWER", "Simes"), 
                      labels=c("Balanced", "Linear"))
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)

kc <- as.character(datC$kMax)
kc[which(datC$kMax==m)] <- "m"
kc[which(datC$kMax==m/2)] <- "m/2"
kc[which(datC$kMax==2*(1-datC$pi0)*m)] <- "2m1"
datC$kMaxC <- kc

datC$ff <- datC$flavor

figName <- sprintf("%s,%s", sname0, simFlavor)
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names
if (stratif) {
    pname2 <- sprintf("%s,SMC", pname2)
}

## plot various statistics vs nominal JFWER
library("ggplot2")

kMaxs <- unique(datC$kMaxC)
families <- unique(datC$family)
confs <- expand.grid(alpha=alphas, kMax=kMaxs, family=families, stringsAsFactors=FALSE)

## - - - - - - - - - - - - - - -
## linear family
## - - - - - - - - - - - - - - -
kk <- "m"
aa <- alphas[1]
fam <- "Linear"

ftag <- sprintf("family=%s,alpha=%s,kMax=%s", fam, 100*aa, kk)

filename <- sprintf("%s,%s,%s.pdf", figName, pname2, ftag)
pathname <- file.path(figPath, filename)
datI <- subset(datC, family==fam & alpha==aa & kMaxC==kk & flavor != "Oracle2")

pal <- RColorBrewer::brewer.pal(4,"YlGnBu")
pal <- RColorBrewer::brewer.pal(5,"BuPu")[-1]
#pal <- RColorBrewer::brewer.pal(5,"OrRd")[-1]
pal <- pal[c(3,2,1,4)]
datI$ff <- plyr::revalue(datI$ff, c("unadjusted"="Simes"))
datI <- subset(datI, pi0 < 0.999)
dat.vline <- expand.grid(rho = 0, pi0 = pi0s[pi0s < 0.999])

p <- ggplot(datI, aes_string(x="SNR", y="value", group="ff", color="ff"))
p <- p + geom_line()
p <- p + facet_grid(pi0 ~ rho,
                    #    scales="free_y",
                    labeller=label_bquote(
                        rows= pi[0]==.(pi0),
                        cols= rho==.(rho)))
#    p <- p + scale_y_continuous(limits = c(0, 0.3))
p <- p + labs(color=sprintf("%s family", fam))
p <- p + geom_point()
p <- p + labs(y="Empirical JER", x=expression(bar(mu)))
p <- p + scale_colour_manual(values = rev(pal))
p <- p + theme(axis.title=element_text(size=16),
               strip.text = element_text(size=12),
               legend.text = element_text(size=12),
               legend.title = element_text(size=12))
p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
p <- p + geom_vline(aes(xintercept=sqrt(-4*log(1-pi0))), linetype="dashed", 
                    data=dat.vline, col="gray")

tbl <- table(datI$pi0, datI$rho)
aspectRatio <- nrow(tbl)/(ncol(tbl) + 0.5)  ## 0.5 gives some slack for legend and y label

pdf(pathname, width=8, height=8*aspectRatio)
print(p)
dev.off()


## - - - - - - - - - - - - - - -
## balanced family
## - - - - - - - - - - - - - - -
kk <- c("m", "10")
aa <- alphas[1]
fam <- "Balanced"

kkt <- paste(kk, collapse="+")
ftag <- sprintf("family=%s,alpha=%s,kMax=%s", fam, 100*aa, kkt)

filename <- sprintf("%s,%s,%s.pdf", figName, pname2, ftag)
pathname <- file.path(figPath, filename)
datI <- subset(datC, family==fam & alpha==aa & kMaxC %in% kk & flavor != "Oracle2")


datI$gg <- sprintf("%s (K=%s)", datI$ff, datI$kMaxC)


pal <- RColorBrewer::brewer.pal(6,"BrBG")
pal <- pal[c(2, 5, 3, 4, 1, 6)]
datI <- subset(datI, flavor != "unadjusted")
datI <- subset(datI, pi0 < 0.999)

p <- ggplot(datI, aes_string(x="SNR", y="value", group="gg", color="gg"))
p <- p + geom_line()
p <- p + facet_grid(pi0 ~ rho,
                    #    scales="free_y",
                    labeller=label_bquote(
                        rows= pi[0]==.(pi0),
                        cols= rho==.(rho)))
#    p <- p + scale_y_continuous(limits = c(0, 0.3))
p <- p + labs(color=sprintf("%s family", fam))
p <- p + geom_point()
p <- p + labs(y="Empirical JER", x=expression(bar(mu)))
p <- p + scale_colour_manual(values = rev(pal))
p <- p + theme(axis.title=element_text(size=16),
               strip.text = element_text(size=12),
               legend.text = element_text(size=12),
               legend.title = element_text(size=12))
p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")

tbl <- table(datI$pi0, datI$rho)
aspectRatio <- nrow(tbl)/(ncol(tbl) + 0.75)

pdf(pathname, width=8, height=8*aspectRatio)
print(p)
dev.off()
