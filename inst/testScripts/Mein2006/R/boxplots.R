boxplots <- function(m, pi0, n, rho, nbSimu, path="resData", figPath="figPath") {
  ## nulls and alternatives
  m0 <- round(m*pi0)
  m1 <- m-m0
  
  simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,nbSimu=%s", m, pi0, n, rho, nbSimu)
  filename <- sprintf("Mein2006,S,%s.xdr", simTags)
  pathname <- file.path(path, filename)
  if (!file.exists(pathname)) {
    print(sprintf("File not found: %s", pathname))
    return()
  }
  S <- loadObject(pathname)
  
  filenameBar <- sprintf("Mein2006,Sbar,%s.xdr", simTags)
  pathnameBar <- file.path(path, filenameBar)
  if (!file.exists(pathname)) {
    print(sprintf("File not found: %s", pathname))
    return()
  }
  Sbar <- loadObject(pathnameBar)

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
