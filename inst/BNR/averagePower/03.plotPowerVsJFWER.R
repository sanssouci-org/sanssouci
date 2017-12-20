filename <- sprintf("%s,%s.rds", sname, pname)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
head(dat)

powerz <-  c("detPow1"="P(S(R,H1)>1", "detPow"="P(S(R,H)>1", 
             "estPow1"="E(S(R,H1))/m1", "estPow"="E(S(R,H))/m1",
             "powBH5"="Power(BH(0.05))", "powBH50"="Power(BH(0.5))", "pow0"="Power({p <= 0.05})",
             "JR"="JFWER")

powerz <-  c("estPow"="(a) R0 = Nm", "powBH5"="(c) R0 = BH(0.05)", "powBH50"="BH(0.5)", "pow0"="(b) R0 = {p <= 0.05}")

gat <- tidyr::gather(dat, "criterion", "value", JR, detPow, detPow1, v0, estPow, estPow1, powBH5, powBH50, pow0)
gat <- subset(gat, criterion %in% names(powerz))

## some reshaping
kc <- as.character(gat$kMax)
kc[which(gat$kMax==m)] <- "m"
kc[which(gat$kMax==m/2)] <- "m/2"
kc[which(gat$kMax==2*(1-gat$pi0)*m)] <- "2m1"
gat$kMaxC <- kc

## some more reshaping
gat$alpha <- as.numeric(gat$alpha)
alphas <- unique(gat$alpha)
alphas <- alphas[which(alphas<=0.25)]

## "balanced" kernel
datB <- subset(gat, flavor == "Step down" & family=="kFWER")
levk <- c("10", "2m1", "m")
datB<- subset(datB, kMaxC %in% levk)
datB$ff <- sprintf("Balanced (K=%s)", datB$kMaxC)

## "linear" kernel
## datS <- subset(gat, family=="Simes" & flavor %in% c("unadjusted", "Single Step", "Step down") & kMax==m)
datS <- subset(gat, family=="Simes" & flavor %in% c("unadjusted", "Step down") & kMax==m)
datS$ff <- factor(datS$flavor, labels=c("Linear", "Simes"))
# datS <- subset(gat, family=="Simes" & flavor %in% c("Step down") & kMax==m)
# datS$ff <- "Linear"

## colors
pal <- RColorBrewer::brewer.pal(7,"PRGn")
pal <- RColorBrewer::brewer.pal(6,"BrBG")
pal <- pal[c(1, 2, 3, 6, 5)]

datC <- rbind(datB, datS)


figName <- sname0
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names

## plot various statistics vs nominal JFWER
library("ggplot2")
x <- "alpha"
rhos <- unique(dat$rho)
rhos <- unique(dat$rho)
sfs <- unique(dat$sf)
confs <- expand.grid(rho=rhos, sf=sfs, x=x, stringsAsFactors=FALSE)

for (ii in 1:nrow(confs)) {
    rr <- confs[ii, "rho"]
    xx <- confs[ii, "x"] 
    ss <- confs[ii, "sf"]
    ftag <- sprintf("rho=%s,SNR=%sx", rr, ss)
    
    
    filename <- sprintf("%s,BalancedVsLinear,%s,%s.pdf", figName, pname2, ftag)
    pathname <- file.path(figPath, filename)
    datI <- subset(datC, rho==rr & sf==ss)
    ## select relevant configs to plot
    datI <- subset(datI,  pi0!=0.999 
                   & criterion!="powBH50"
                   & (family!="Linear" | kMax==m))
    datI$plab <- sprintf("pi[0]=%s ; mu=%s", datI$pi0, round(getSNR(datI$pi0, ss), 2))
    datI$plab <- "pi[0]"
    p <- ggplot(datI, aes_string(x="alpha", y="value", group="ff", color="ff"))
    p <- p + geom_line()
    p <-
        p + facet_grid(criterion ~ pi0+SNR,
                        scales="free_y",
                        labeller=label_bquote(rows=.(powerz[[criterion]]),
                                              cols= pi[0]:.(pi0)~-~bar(mu):.(round(SNR, 1))))
    p <- p + scale_x_continuous(breaks=round(alphas, 2), minor_breaks=NULL, limits = range(alphas))
    p <- p + scale_y_continuous(minor_breaks=NULL, limits=c(0,1))
    ##    p <- p + scale_y_continuous(minor_breaks=NULL)
    p <- p + theme(axis.text.x=element_text(angle=90))
    p <- p + labs(color="Family",
                  linetype=expression(lambda-adjustment))
    p <- p + geom_point()
    p <- p + labs(y="Averaged power", x="Target JER level")
    p <- p + scale_colour_manual(values = pal)
    p <- p + theme(axis.title=element_text(size=16),
                   strip.text = element_text(size=12),
                   legend.text = element_text(size=12),
                   legend.title = element_text(size=12))
    pdf(pathname, width=8, height=7)
    print(p)
    dev.off()
}