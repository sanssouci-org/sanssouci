library("R.utils")
library("Hmisc")

study <- "Mein2006,sansSouci,0.0.5"
path <- file.path("resData", study)
path <- Arguments$getReadablePath(path)

alpha <- 0.05
n <- 100

tag <- sprintf("n=%s,alpha=%s", n, alpha)

filename <- sprintf("resM,%s.xdr", tag)
pathname <- file.path(path, filename)
mea <- loadObject(pathname)

filename <- sprintf("resS,%s.xdr", tag)
pathname <- file.path(path, filename)
se <- loadObject(pathname)

x <- as.numeric(colnames(mea))
ylim <- c(min(mea-2*se), max(mea+2*se))
ylim[1] <- 0
ylim[2] <- alpha*1.2
xlim <- range(x)
nr <- nrow(mea)
cols <- 1:nr
pchs <- 15:17
rn <- as.numeric(rownames(mea))
lgd <- sapply(1:length(rn), function(ii) {
  as.expression(substitute(rho*"="*rr, list(rr=rn[ii])))
})

figPath <- file.path("fig", study)
figPath <- Arguments$getWritablePath(figPath)
figName <- sprintf("Mein2006,%s", tag)
figName <- gsub("\\.", "_", figName)
filename <- sprintf("%s.pdf", figName)
pathname <- file.path(figPath, filename)

pdf(pathname)
plot(NA, ylim=ylim, xlim=xlim, xlab=expression(pi[0]), ylab=expression(hat(p)))
for (ii in 1:nr) {
  y <- mea[ii, ]
  yp <- y+2*se[ii, ]
  ym <- y-2*se[ii, ]
  errbar(x, y, yp, ym, col=ii, errbar.col=ii, add=TRUE, t='b', pch=pchs[ii])
}
abline(h=alpha, col="gray", lty=2)
legend("bottomright", legend=lgd, col=cols, lty=1, pch=pchs, horiz=FALSE, box.lwd=0)
dev.off()
