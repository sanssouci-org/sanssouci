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
##dat$estPowM1 <- dat$estPow/(1-dat$pi0)

powerz <-  list("Power1"="P(S(R,H1)>1", "Power"="P(S(R,H)>1", "s1"="E(S(R,H1))/m1", "s1"="E(S(R,H1))/m1")
##                res <- rbind(JR=rej0>0, detPow1=rej1>0, detPow=rej01>0, v0, estPow1=s1, estPow=s01)
powerz <-  list("detPow1"="P(S(R,H1)>1", "detPow"="P(S(R,H)>1", 
                "estPow1"="E(S(R,H1))/m1", "estPow"="E(S(R,H))/m1",
                "powBH"="half of BH(0.05)", "pow0"="half of {p <= 0.05}")

## some reshaping
#datC <- subset(dat, kMax %in% as.character(c(10, 100, m, m/2))  & flavor != "unadjusted")
datC <- subset(dat, flavor != "unadjusted")

mm <- grep("Oracle|alpha", datC$flavor)
datC <- datC[-mm, ]
levels(datC$family) <- list("Balanced"="kFWER", "Linear"="Simes")
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)
datC$ff <- factor(paste(datC$family, datC$flavor))
datC$ff <- factor(paste(datC$family, datC$kMaxC))

if (FALSE) {
    ## drop outliers?
    ww <- which(with(datC, dep==0 & pi0==0.8 & alpha==0.025 & Power < 0.5))
    datC <- datC[-ww, ]
}

figName <- sname0
figPath <- file.path("fig", sprintf("BNR,%s", ptag))
figPath <- R.utils::Arguments$getWritablePath(figPath)
pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names

## plot Power vs JFWER
library("ggplot2")

x <- c("JR", "alpha")[2]
SNRs <- unique(dat$SNR)
#SNRs <- 2

confs <- expand.grid(SNR=SNRs, power=names(powerz), x=x, stringsAsFactors=FALSE)

for (ii in 1:nrow(confs)) {
    snr <- confs[ii, "SNR"]
    pp <- confs[ii, "power"] 
    xx <- confs[ii, "x"] 
    pow <- powerz[[pp]]
    ftag <- sprintf("SNR=%s,power=%s", snr, gsub("/", "--", pow))
    if (xx=="JR") {
        filename <- sprintf("PowerVsJFWER,%s,%s,%s.pdf", figName, pname2, ftag)
    } else if (xx=="alpha") {
        filename <- sprintf("PowerVsNominalJFWER,%s,%s,%s.pdf", figName, pname2, ftag)
    } else {
        stop("I don't know what to do when x=", xx)
    }
        
    pathname <- file.path(figPath, filename)

    datI <- subset(datC, SNR==snr)
    
    pdf(pathname)
    p <- ggplot(datI, aes_string(x=xx, y=pp, group="ff", color="ff"))
    ##p <- p + geom_line() + geom_point(aes(shape=factor(alpha))) + scale_shape_identity()
    #p <- p + geom_line(aes(linetype=flavor))
    p <- p + geom_line()
    ##p <- p + facet_grid(pi0 ~ rho, scales="free_y")
    p <- p + facet_grid(pi0 ~ dep,
                        scales="free_y",
                        labeller=label_bquote(
                            cols= rho==.(dep),
                            rows= pi[0]==.(pi0)))
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
    p <- p + labs(y=pow)
    p <- p + scale_color_brewer(type="div")
    print(p)
    dev.off()
}

## sanity check wrt power definitions
pow1 <- subset(dat, family=="kFWER" & flavor =="Single Step")
pow1sd <- subset(dat, family=="kFWER" & flavor =="Step down")

identical(pow1$Power, pow1sd$Power)

identical(pow1$Power1, pow1sd$Power1)  
table(sign(pow1sd$Power1-pow1$Power1))

identical(pow1$s1, pow1sd$s1)
table(sign(pow1sd$s1-pow1$s1))
