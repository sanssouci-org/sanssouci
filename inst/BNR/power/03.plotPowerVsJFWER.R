library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")
plan(remote, workers = rep("bernoulli", 10))

## retrieve results (remote!)
dat %<-% {
    rpath <- file.path("~/Documents/Packages/sanssouci", path)
    fls <- list.files(rpath, full.names=TRUE)
    id <- gsub(".rds$", "", basename(fls))
    names(fls) <- id
    plyr::ldply(fls, readRDS, .id="id")
}
head(dat)

## some reshaping
datC <- subset(dat, kMax==as.character(m) & SNR==2 & flavor != "unadjusted")
mm <- grep("Oracle|alpha", datC$flavor)
datC <- datC[-mm, ]
levels(datC$family) <- list("Balanced"="kFWER", "Simes"="Simes")
datC$alpha <- as.numeric(datC$alpha)
alphas <- unique(datC$alpha)
datC$ff <- factor(paste(datC$family, datC$flavor))

if (FALSE) {
    ## drop outliers?
    ww <- which(with(datC, dep==0 & pi0==0.8 & alpha==0.025 & Power < 0.5))
    datC <- datC[-ww, ]
}

pname2 <- gsub("0\\.", "", pname)  ## to avoid '.' in LaTeX file names
filename <- sprintf("PowerVsJFWER,%s,%s.pdf", figName, pname2)
pathname <- file.path(figPath, filename)

## plot Power vs JFWER
## caution: we're interested in 's1' (ie average power), not 'Power' (ie detection power)

pdf(pathname)
p <- ggplot(datC, aes(x=JR, y=s1, group=ff, color=family))
##p <- p + geom_line() + geom_point(aes(shape=factor(alpha))) + scale_shape_identity()
p <- p + geom_line(aes(linetype=flavor))
p <- p + geom_point(aes(shape=factor(alpha))) + scale_shape_manual(values=1:9)
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
              linetype=expression(lambda-adjustment),
              shape="Target JFWER level")
print(p)
dev.off()
