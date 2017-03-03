filename <- sprintf("%s,%s.rds", sname, pname)
pathname <- file.path(resPath, filename)
dat <- readRDS(pathname)
head(dat)

gat <- tidyr::gather(dat, "criterion", "value", JR, detPow, detPow1, v0, estPow, estPow1, powBH5, powBH50, pow0)

## some reshaping
datC <- subset(gat, flavor == "Step down" & criterion=="detPow" & beta !=2/3)
datC$family <- factor(datC$family, levels=c("kFWER", "Simes"), labels=c("Balanced", "Linear"))
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)
alphas <- alphas[which(alphas<=0.25)]

kc <- as.character(datC$kMax)
kc[which(datC$kMax==m)] <- "m"
kc[which(datC$kMax==m/2)] <- "m/2"
kc[which(datC$kMax==2*(1-datC$pi0)*m)] <- "2m1"
datC$kMaxC <- kc

levk <- c("10", "m")
datC <- subset(datC, kMaxC %in% levk & alpha %in% alphas)
levs <- c(paste("Balanced", levk), paste("Linear", rev(levk)))
#datC$ff <- factor(paste(datC$family, datC$kMaxC), levels=levs)
datC$ff <- sprintf("%s (K=%s)", datC$family, datC$kMaxC)
#datC$ff <- factor(ff, levels=levs)

figName <- sname0
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names

## plot various statistics vs nominal JFWER
library("ggplot2")

x <- "alpha"
rhos <- unique(dat$rho)
confs <- expand.grid(rho=rhos, x=x, stringsAsFactors=FALSE)

for (ii in 1:nrow(confs)) {
    rr <- confs[ii, "rho"]
    xx <- confs[ii, "x"] 
    ftag <- sprintf("rho=%s", rr)
    
    filename <- sprintf("%s,BalancedVsLinear,%s,%s.pdf", figName, pname2, ftag)
    pathname <- file.path(figPath, filename)
    datI <- subset(datC, rho==rr)
#    datI$rf <- factor(datI$r, levels=sort(unique(datI$r), decr=TRUE))
#    datI$rf <- factor(datI$r, levels=sort(unique(datI$r), decr=FALSE))
    
    p <- ggplot(datI, aes_string(x=xx, y="value", group="ff", color="ff"))
    p <- p + geom_line()
    p <- p + facet_grid(r ~ beta,
                        scales="free_y",
                        labeller=label_bquote(
                            rows= .(r),
                            cols= beta==.(round(beta, 2))),
                        as.table=FALSE)
    p <- p + scale_x_continuous(breaks=round(alphas, 2), minor_breaks=NULL, limits = range(alphas))
    p <- p + scale_y_continuous(minor_breaks=NULL, limits=c(0,1))
    ##    p <- p + scale_y_continuous(minor_breaks=NULL)
    p <- p + theme(axis.text.x=element_text(angle=90))
    p <- p + labs(color="Family",
                  linetype=expression(lambda-adjustment))
    p <- p + geom_point()
    p <- p + labs(y="Detection power", x="Target JER level")
    p <- p + scale_color_brewer(type="div")
    p <- p + theme(axis.title=element_text(size=16),
                   strip.text = element_text(size=12),
                   legend.text = element_text(size=12),
                   legend.title = element_text(size=12))
    pdf(pathname, width=9)
    print(p)
    dev.off()
}

detMax <- function(beta) (1-sqrt(1-beta))^2
detHC <- function(beta) {
    r <- detMax(beta)
    r[beta<3/4] <- beta[beta<3/4]-1/2
    r
}
bs <- unique(datC$beta)
xy <- expand.grid(beta=bs, r=rs)
isDetMax <- (xy$r>=detMax(xy$beta))
isDetHC <- (xy$r>=detHC(xy$beta))
xy$Region <- "Detectable"
xy$Region[!isDetMax & isDetHC] <- "Detectable by HC only"
xy$Region[!isDetHC] <- "Undetectable"

pal <- RColorBrewer::brewer.pal(3,"Dark2")
sz <- 14

q <- qplot(beta, r, data=xy, color=Region, shape=Region)
q <- q + stat_function(fun=detMax, colour = "black")
q <- q + stat_function(fun=detHC, colour = "black")
q <- q + xlab(expression(beta)) + geom_point(size=3)
q <- q +  theme_bw() + theme(panel.grid.major = element_line(colour = "grey90"),
                             panel.grid.minor = element_blank(),
                             #legend.position = "top", legend.direction="vertical", 
                             legend.text=element_text(size=sz),
                             legend.title=element_text(size=sz),
                             axis.title=element_text(size=sz),
                             axis.text=element_text(size=sz),
                             strip.text = element_text(size=sz))
q <- q + scale_x_continuous(breaks=round(bs, 2), minor_breaks=NULL, limits = range(bs))
q <- q + scale_y_continuous(breaks=round(rs, 2), minor_breaks=NULL, limits = range(rs))
q <- q + scale_colour_manual(values = pal)
q <- q + annotate("text", x = 0.555, y = 0.14, label = "FDR/FWER", angle=11)
q <- q + annotate("text", x = 0.555, y = 0.025, label = "HC", angle=22)

pdf("fig/DJ2004.pdf", width=8, height=4)
print(q)
dev.off()
