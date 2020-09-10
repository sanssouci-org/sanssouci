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

