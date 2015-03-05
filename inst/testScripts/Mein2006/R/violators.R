violators <- function(m=1e3, pi0=1, n=100, rho=0, B=1e3, nbSimu=1e3, alpha=0.05, SNR=NA, path="resData") {
  ## nulls and alternatives
  m0 <- round(m*pi0)
  m1 <- m-m0

  if (is.na(SNR)) {
    simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,B=%s,nbSimu=%s,alpha=%s", m, pi0, n, rho, B, nbSimu, alpha)
  } else {
    simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,SNR=%s,B=%s,nbSimu=%s,alpha=%s", m, pi0, n, rho, SNR, B, nbSimu, alpha)
  }    
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
  mvio <- NA
  ## "violators"
  vio <- apply(Sbar-S, 2, FUN=function(x) sum(x>0))
  mvio <- mean(vio>0)
  svio <- sd(vio>0)

  if (FALSE) { ## sanity check
    V <- matList[["V"]]
    Vbar <- matList[["Vbar"]]
    vio2 <- apply(Vbar-V, 2, FUN=function(x) sum(x<0))
    stopifnot(identical(vio, vio2))
  }
  
  list(mean=mvio, sd=svio)
}
