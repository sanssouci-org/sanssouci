boxplotsV <- function(m, pi0, n, rho, nbSimu, path="resData", figPath="figPath", plot=TRUE) {
  ## nulls and alternatives
  m0 <- round(m*pi0)
  m1 <- m-m0
  
  simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,nbSimu=%s", m, pi0, n, rho, nbSimu)
  fTags <- c("S", "Sbar", "V", "Vbar")
  filenames <- sprintf("Mein2006,%s,%s.xdr", fTags, simTags)
  pathnames <- file.path(path, filenames);

  for (pathname in pathnames) {
    if (!file.exists(pathname)) {
      print(sprintf("File not found: %s", pathname))
      return()
    }
  }
  matList <- lapply(pathnames, loadObject)
  names(matList) <- fTags
  S <- matList[["S"]]
  Sbar <- matList[["Sbar"]]
  V <- matList[["V"]]
  Vbar <- matList[["Vbar"]]
  mvio <- NA
  dropped1 <- 0
  if (m1>0) {
    ## need to take the max value of 'Sbar' acroos all identical values in 'S'
    uS <- seq_len(m1)  ## unique values in 'S';
    matSbar <- matrix(NA, nbSimu, m1)
    colnames(matSbar) <- uS
    for (kk in 1:nbSimu) {
      test <- (length(unique(S[, kk]))==m1) ## sanity check
      if (test) {
        ## SbarKK <- aggregate(Sbar[, kk], by=list(S[, kk]), max)$x        
        SbarKK <- Sbar[which(diff(c(S[, kk], m1+1))>=1), kk]  ## heaps faster
        matSbar[kk, ] <- SbarKK
      } else {
        warning("Dropping a simulation with missing values for 'S'")
        dropped1 <- dropped1+1
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

    if (plot) {
      figname <- sprintf("boxplot,S,%s.pdf", simTags)
      pathname <- file.path(figPath, figname)
      pdf(pathname)
      boxplot(matSbar[, idxs], at=idxs, outline=TRUE, axes=FALSE, xlim=xlim, ylim=ylim, boxwex=bwx)
      axis(2, cex.axis=2, las=2)
      axis(1, at=idxs, labels=uS[idxs], cex.axis=1.5, las=2)
      abline(a=0, b=ylim[2]/xlim[2])
      dev.off()
    }

    ## "violators"
    vio <- apply(Sbar-S, 2, FUN=function(x) sum(x>0))
  }
  dropped0 <- 0
  if (m0>0) {
    ## need to take the *min* value of 'Vbar' acroos all identical values in 'V'
    uV <- seq_len(m0+1)  ## unique values in 'V';
    matVbar <- matrix(NA, nbSimu, m0+1);
    colnames(matVbar) <- uV
    for (kk in 1:nbSimu) {
      test <- (length(unique(V[, kk]))==m0+1) ## sanity check
      if (test) {
        ## VbarKK <- aggregate(Vbar[, kk], by=list(V[, kk]), min)$x
        ## VbarKK <- tapply(Vbar[, kk], V[, kk], min)  ## faster (yet not perfect)
        VbarKK <- Vbar[which(diff(c(-1, V[, kk]))>=1), kk]
        matVbar[kk, ] <- VbarKK
      } else {
        warning("Dropping a simulation with missing values for 'V'")
        dropped0 <- dropped0+1
      }
    }

    ## VbarKK <- aggregate(Vbar[, kk], by=list(V[, kk]), min)$x
    ## VbarKK <- tapply(Vbar[, kk], V[, kk], min)
    ## Vbar[which(diff(c(-1, Vkk)>=1), kk]
    
    len <- 20
    if (m0>50) {
      idxs <- round(seq(from=0, to=m0, by=m0/len))
      idxs[1] <- 1
    } else {
      idxs <- 1:m0
    }
    ylim <- c(1, m0)
    xlim <- c(1, m0)
    bwx <- round(m0/length(idxs))*0.8

    if (plot) {
      figname <- sprintf("boxplot,V,%s.pdf", simTags)
      pathname <- file.path(figPath, figname)
      pdf(pathname)
      boxplot(matVbar[, idxs], at=idxs, outline=TRUE, axes=FALSE, xlim=xlim, ylim=ylim, boxwex=bwx)
      axis(2, cex.axis=1.5, las=2)
      axis(1, at=idxs, labels=colnames(Vbar)[idxs], cex.axis=1.4, las=2)
      abline(a=0, b=ylim[2]/xlim[2])
      dev.off()

      toPlot <- c(1,5, 10, 20, 50)
      cols <- ltys <- seq(along=toPlot)
      lgd <- sprintf("%02d/%s", toPlot, m)
      
      o <- order(matVbar[, round(m/2)])
      trajs <-  t(matVbar[o[toPlot], ])
      figname <- sprintf("traj,V,%s.pdf", simTags)
      pathname <- file.path(figPath, figname)
      pdf(pathname)
      boxplot(matVbar[, idxs], at=idxs, outline=TRUE, axes=FALSE, xlim=xlim, ylim=ylim, boxwex=bwx, border="lightgray")
      ## plot(NA, axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="")
      matlines(trajs, t='l', col=cols, lty=ltys)
      ## matplot(trajs, t='l', col=cols, lty=ltys, add=TRUE)
      axis(2, cex.axis=1.5, las=2)
      axis(1, at=idxs, labels=colnames(Vbar)[idxs], cex.axis=1.4, las=2)
      abline(a=0, b=ylim[2]/xlim[2])
      legend("bottomright", lgd, col=cols, lty=ltys, title="Rank of the curve")
      dev.off()
    }
    
    ## "violators"
    vio2 <- apply(Sbar-S, 2, FUN=function(x) sum(x>0))
    if (exists("vio")) {
      stopifnot(identical(vio, vio2))
    } else {
      vio <- vio2
    }
  }
  mvio <- mean(vio>0)
  print(mvio)
  svio <- sd(vio>0)
  list(mean=mvio, sd=svio, dropped0=dropped0, dropped1=dropped1)
}
