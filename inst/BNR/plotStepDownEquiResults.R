path <- "../../../../resData/stepDownEqui,Rcpp"
typeOfSNR <- c("constantSNR", "Pareto")[1]

res <- NULL
fls <- list.files(path)
str(fls)
for (ff in fls) {
    pathname <- file.path(path, ff)
    dat <- readRDS(pathname)
    resMat <- Reduce(rbind, dat)    
    res <- rbind(res, colMeans(resMat>0))
}
colnames(res) <- c("Single step", "Step down", "Oracle")
## lgd <- gsub("stepDownEqui,m=1000,(.*),B.*", "\\1", fls)
lgd <- gsub("stepDownEqui,m=1000,(.*).rds", "\\1", fls)
## lgd <- gsub("stepDownEqui,m=1000,(.*),B=10000,.*,nbSimu=1000.rds", "\\1", fls)
rownames(res) <- lgd
print(res)

pattern <- "pi0=(.*),rho=(.*),SNR=(.*),B=(.*),alpha=(.*),nbSimu=(.*)"
pi0 <- gsub(pattern, "\\1", lgd)
pi0 <- as.numeric(gsub("_", "\\.", pi0))

rho <- gsub(pattern, "\\2", lgd)
rho <- as.numeric(gsub("_", "\\.", rho))

SNRc <- gsub(pattern, "\\3", lgd)
SNR <- as.numeric(SNRc)

if (typeOfSNR=="constantSNR") {
    ww <- which(!is.na(SNR))
} else if (typeOfSNR=="Pareto") {
    ww <- which(is.na(SNR))
    SNR <- SNRc
}

B <- gsub(pattern, "\\4", lgd)
B <- as.numeric(B)

alpha <- gsub(pattern, "\\5", lgd)
alpha <- as.numeric(gsub("_", "\\.", alpha))

nbSimu <- gsub(pattern, "\\6", lgd)
nbSimu <- as.numeric(nbSimu)

colnames(res) <- c("SingleStep", "StepDown", "Oracle")
dat <- cbind(data.frame(pi0=pi0, rho=rho, SNR=SNR, B=B, alpha=alpha, nbSimu=nbSimu, stringsAsFactors=FALSE), res)
rownames(dat) <- NULL
dat <- dat[ww, ]  ## keep either Pareto or non-Pareto results

#datS <- dat
datS <- subset(dat, nbSimu==10000 & SNR != 10)
#datS <- subset(dat, nbSimu==2000)
#datS <- subset(dat, B==10000 & nbSimu==2000)
#datS <- subset(dat, B==5000 & nbSimu==1000)

alpha <- unique(datS$alpha)
n <- unique(datS$nbSimu)

if (length(n) != 1L) {
    warning("Mixing several simulation settings here!")
    n <- n[1]
}

datSC <- datS[, -match(c("B", "alpha", "nbSimu"), names(datS))]

## table
library("xtable")
tab <- xtable(datSC)
o <- order(datSC$rho, datSC$pi0, datSC$SNR)
print(tab[o, ], include.rownames=FALSE)
datSC$SNR

library("R.utils")
figPath <- "../../../../fig/BNR"
figPath <- Arguments$getWritablePath(figPath)

## plot
library("ggplot2")
p <- ggplot(datSC, aes(x=SNR, y=SingleStep, group=interaction(rho, pi0), color=rho, lty=as.factor(pi0)))
p <- p + geom_line() + geom_point()

library("reshape2")
d <- melt(datSC, id.vars=c("pi0", "rho", "SNR"), value.name="Empirical JFWER achieved", variable.name="Method")
d$se <- sqrt(d$`Empirical JFWER achieved` * (1-d$`Empirical JFWER achieved`))/sqrt(n)
limits <- aes(ymax=`Empirical JFWER achieved` + se, ymin=`Empirical JFWER achieved`-se)

p <- ggplot(d, aes(x=SNR, y=`Empirical JFWER achieved`, group=Method, color=Method))
p <- p + geom_line() + geom_point() + geom_errorbar(limits, width=0.1)
p <- p + facet_grid(pi0 ~ rho)
#p <- p + facet_grid(pi0 ~ rho, labeller=labeller(rho=label_bquote(rho*"="*.(x)), pi0=label_bquote(pi[0]*"="*.(x))))
if (typeOfSNR=="constantSNR") {
    p <- p + scale_x_continuous(breaks=unique(d$SNR)) + xlab(expression(mu))
} else {
    p <- p + scale_x_discrete(breaks=unique(d$SNR)) + xlab(expression(mu))
    p <- p + theme(axis.text.x=element_text(angle=90))
}
p <- p + scale_y_continuous(breaks=c(0.1, 0.2, 0.25))
p <- p + geom_hline(aes(yintercept=alpha), linetype="dashed")
p <- p + geom_hline(aes(yintercept=alpha*pi0), linetype="dotted")
p


filename <- "stepDownEqui.pdf"
if (typeOfSNR=="Pareto") {
    filename <- "stepDownEqui,Pareto.pdf"
}
pathname <- file.path(figPath, filename)
ggsave(pathname)

# bof:
if (typeOfSNR == "Pareto") {
    library("actuar")
    filename <- "Pareto.pdf"
    ppar <- unique(datS$SNR)
    pattern <- "Pareto\\(([0-9]),([0-9]),([0-9]))"
    xmin <- as.numeric(gsub(pattern, "\\1", ppar))
    shape <- as.numeric(gsub(pattern, "\\2", ppar))
    scale <- as.numeric(gsub(pattern, "\\3", ppar))
    q <- qplot(c(0.01, 10), geom = 'blank')
    mat <- NULL
    pathname <- file.path(figPath, filename)
    pdf(pathname)
    for (pp in seq(along = ppar)) {
        curve(dpareto(x + xmin[pp], shape = shape[pp], scale = scale[pp]), add=pp>1, from=0, to=100,  col=pp, log="y")
    }
    dev.off()
}
