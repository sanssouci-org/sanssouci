## Retrieving results
flavor <- c("2006", "2013")[1]
study <- sprintf("Mein%s", flavor)

path <- file.path("resData", study)
path <- Arguments$getReadablePath(path)

filenames <- list.files(path)
pathnames <- file.path(path, filenames)

figPath <- file.path("fig", study, "fig1")
figPath <- Arguments$getWritablePath(figPath)

boxplots <- function(m, pi0, n, rho, nbSimu) {
  simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,nbSimu=%s", m, pi0, n, rho, nbSimu)

  ## nulls and alternatives
  m0 <- round(m*pi0)
  m1 <- m-m0
  
  filename <- sprintf("Mein2006,S,%s.xdr", simTags)
  pathname <- file.path(path, filename)
  if (!file.exists(pathname)) {
    return()
  }
  S <- loadObject(pathname)
  
  filenameBar <- sprintf("Mein2006,Sbar,%s.xdr", simTags)
  pathnameBar <- file.path(path, filenameBar)
  if (!file.exists(pathname)) {
    return()
  }
  Sbar <- loadObject(pathnameBar)

  if (FALSE) {
    plot(S[, 1], Sbar[, 1], t='l', xlab=expression(S(t)), ylab=expression(underline(S)(t)))
    abline(a=0, b=1, col=2)
    points(S[, 1], Sbar[, 1], t='l', xlab=expression(S(t)), ylab=expression(underline(S)(t)))
    ## apply(S, 2, FUN=function(x) length(unique(x)))
  }
  
  uS <- seq_len(m1)  ## unique values in 'S';
  matSbar <- matrix(NA, nbSimu, m1)
  colnames(matSbar) <- uS
  for (kk in 1:nbSimu) {
    test <- (length(unique(S[, kk]))==m1) ## sanity check
    if (test) {
      SbarKK <- aggregate(Sbar[, kk], by=list(S[, kk]), max)$x
      matSbar[kk, ] <- SbarKK
    } else {
      warning("")
    }
  }

  len <- 20
  if (m1>50) {
    idxs <- round(seq(from=0, to=m1, by=m1/len))
    idxs[1] <- 1
  } else {
    idxs <- 1:m1
  }
  ylim <- c(1, m1)
  xlim <- c(1, m1)
  bwx <- round(m1/length(idxs))*0.8
  
  figname <- sprintf("boxplot,%s.pdf", simTags)
  pathname <- file.path(figPath, figname)
  pdf(pathname)
  boxplot(matSbar[, idxs], at=idxs, outline=TRUE, axes=FALSE, xlim=xlim, ylim=ylim, boxwex=bwx)
  axis(2, cex.axis=2, las=2)
  axis(1, at=idxs, labels=uS[idxs], cex.axis=1.5, las=2)
  abline(a=0, b=ylim[2]/xlim[2])
  dev.off()

  ## "violators"
  vio <- apply(Sbar-S, 2, FUN=function(x) sum(x>0))
  print(mean(vio>0))  
}

nbSimu <- 1000
ms <- c(100, 1000)[2]
rhos <- c(0, 0.4)
ns <- c(20, 40, 60, 80)
pi0s <- c(0.7, 0.6)[1]
for (m in ms) for (pi0 in pi0s) for (n in ns) for (rho in rhos) boxplots(m, pi0, n, rho, nbSimu)
