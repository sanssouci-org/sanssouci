## For Figures 2 and 3 in the BNR paper
library("future")
library("listenv")
##computeNodes <- c("cauchy", "leibniz", "bolzano", "shannon", "euler", "hamming", "bernoulli")

library("sansSouci")

## setup
m <- 1e3
kMaxs <- c(1, 2, 10, m)
B <- 1e4
nbSimu <- 100

alpha <- 0.2

rhos <- seq(from=0, to=1, by=0.1)

oneSim <- function(m, rho, B, alpha, pi0=1, SNR=0, kMaxs=kMaxs) {
    sim <- sansSouci::simulateEqui(m, rho, B, pi0=pi0, SNR=SNR)
    configs <- expand.grid(kMax=kMaxs, refFam=c("kFWER", "Simes"), stringsAsFactors=FALSE)
    lambdas <- numeric(0L)
    for (rr in 1:nrow(configs)) {
        fam <- configs[rr, "refFam"]
        kMax <- configs[rr, "kMax"]
        resJ <- jointFWERControl(sim$X0, refFamily=fam, alpha, kMax=kMax, verbose=FALSE)
        lambdas <- c(lambdas, resJ$lambda)
    }
    df <- cbind(configs, lambda=lambdas)
    return(df)
}

#plan(remote, workers = rep("bernoulli", 100))
plan(multiprocess, workers = 150)
options("future.wait.times"=1e4)

resr <- listenv()
for (rr in seq(along=rhos)) {
    rho %<-% rhos[rr]
    res %<-% listenv()
    for (ss in 1:nbSimu) {
        if (ss %% 20==0) { print(ss);}
        res[[ss]] %<-% {
            library("sansSouci")
            oneSim(m, rho, B, alpha, kMaxs=kMaxs)
        }
    }
    names(res) %<-% 1:nbSimu
    dat %<-% plyr::ldply(as.list(res), data.frame, .id="sid")
    resr[[rr]] %<-% dat
}
names(resr) <- rhos
lambdas <- ldply(resr, data.frame, .id="rho")

#saveRDS(lambdas, file="resData/lambdas.rds")
#lambdas <- readRDS("resData/lambdas.rds")

ptag <- sprintf("sansSouci_%s", packageVersion("sansSouci"))
date <- Sys.Date()
figPath <- file.path("fig", sprintf("BNR,%s,%s", ptag, date))
figPath <- R.utils::Arguments$getWritablePath(figPath)
figName <- "lambda"
library("ggplot2")

for (fam in c("kFWER", "Simes")) {
    filename <- sprintf("%s,%s.pdf", figName, fam)
    pathname <- file.path(figPath, filename)
    
    dat <- subset(lambdas, refFam==fam)
    dat$lambda <- dat$lambda/alpha
    mns <- aggregate(dat$lambda, list(rho=dat$rho, kMax=dat$kMax), mean)
    names(mns)[3] <- "lambda"
    mns$K <- factor(mns$kMax)
    
    pdf(pathname, width=6, height=3)
    p <- ggplot(mns, aes(rho, lambda, color=K, group=K))
    p <- p + xlab(expression(rho))
    p <- p + ylab(expression(lambda(alpha, N[m])/alpha))
    p <- p + geom_point()
    p <- p + geom_line()    
    if (fam=="Simes") {
        p <- p + coord_cartesian(ylim=c(1, 10))
    }
    p <- p + theme_bw() + theme(panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_blank()
        )
    print(p)
    dev.off()
}
