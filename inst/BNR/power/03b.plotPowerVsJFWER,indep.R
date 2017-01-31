library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(remote, workers = rep("bernoulli", 2)) ## retrieve results (remote!)
plan(eager) ## retrieve results (local!)

head(dat)

gat <- tidyr::gather(dat, "criterion", "value", JR, detPow, estPow, powBH5, powBH50, pow0)

powerz <-  list("detPow1"="P(S(R,H1)>1", "detPow"="P(S(R,H)>1", 
                "estPow1"="E(S(R,H1))/m1", "estPow"="E(S(R,H))/m1",
                "powBH5"="half of BH(0.05)", "powBH50"="half of BH(0.5)", "pow0"="half of {p <= 0.05}",
                "JR"="JFWER")

## some reshaping
#datC <- subset(dat, kMax %in% as.character(c(10, 100, m, m/2))  & flavor != "unadjusted")
datC <- subset(gat, flavor == "Step down")
datC$family <- factor(datC$family, levels=c("kFWER", "Simes"), labels=c("Balanced", "Linear"))
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)

levk <- as.character(sort(as.numeric(unique(datC$kMax))))
levs <- c(paste("Balanced", levk), paste("Linear", rev(levk)))
datC$ff <- factor(paste(datC$family, datC$kMax), levels=levs)

figName <- sname0
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names

## plot various statistics vs nominal JFWER
library("ggplot2")

x <- "alpha"
SNRs <- unique(dat$SNR)
#SNRs <- 2

confs <- expand.grid(SNR=SNRs, x=x, stringsAsFactors=FALSE)

for (ii in 1:nrow(confs)) {
    snr <- confs[ii, "SNR"]
    xx <- confs[ii, "x"] 
    ftag <- sprintf("SNR=%s", snr)
    if (xx=="JR") {
        filename <- sprintf("BalancedVsLinear,indep,%s,%s,%s.pdf", figName, pname2, ftag)
    } else if (xx=="alpha") {
        filename <- sprintf("BalancedVsLinear,indep,%s,%s,%s.pdf", figName, pname2, ftag)
    } else {
        stop("I don't know what to do when x=", xx)
    }
        
    pathname <- file.path(figPath, filename)
    datI <- subset(datC, SNR==snr & kMax %in% c(10, 20, m))
    
    pdf(pathname)
    p <- ggplot(datI, aes_string(x=xx, y="value", group="ff", color="ff"))
    p <- p + geom_line()
    p <- p + facet_grid(criterion ~ pi0,
                        scales="free_y",
                        labeller=label_bquote(
                            rows= .(powerz[[criterion]]),
                            cols= pi[0]==.(pi0)))
    vdat <- data.frame(x=unique(datC$alpha))
    #p <- p + geom_vline(aes(xintercept=x), data=vdat, linetype="dashed", color="gray")
    p <- p + scale_x_continuous(breaks=round(alphas, 2), minor_breaks=NULL, limits = range(alphas))
    p <- p + scale_y_continuous(minor_breaks=NULL)
    p <- p + theme(axis.text.x=element_text(angle=90))
    p <- p + labs(color="Family",
                  linetype=expression(lambda-adjustment))
    if (xx=="JR") {  ## not useful otherwise
        p <- p + geom_point(aes(shape=factor(alpha))) + scale_shape_manual(values=1:9)
        p <- p + labs(shape="Target JR level")
    } else {
        p <- p + geom_point()
    }
    p <- p + labs(y="criterion")
    p <- p + scale_color_brewer(type="div")
    print(p)
    dev.off()
}

