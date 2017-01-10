runMeinshausen2006 <- function(m, pi0, n, rho, nbSimu, p=0.5, B=1e4, alpha=0.05, SNR=1,
                               flavorH=c("2006", "2013", "sansSouci"), ## flavor for the "howmany" package
                               sort=FALSE,  ## sort=TRUE fixes a bug in 'howmany'...
                               path="resData",
                               overwrite=FALSE,
                               mc.cores=2L, ## number of cores to be used by 'parallel::mclapply'
                               verbose=FALSE) {
  if (verbose) {
    print(flavorH)
    print(sprintf("n=%s", n))
    print(sprintf("rho=%s", rho))
    print(sprintf("B=%s", B))
    print(sprintf("nbSimu=%s", nbSimu))
    step <- max(round(nbSimu/100), 1)
  }
  ## nulls and alternatives
  m0 <- round(m*pi0)
  m1 <- m-m0
  H0 <- 1:m0
  if (m0 < m) {
    H1 <- seq(from=m0+1, to=m, by=1)
  } else {
    H1 <- integer(0)
  }
  
  simTags <- sprintf("m=%s,pi0=%s,n=%s,rho=%s,SNR=%s,B=%s,nbSimu=%s,alpha=%s", m, pi0, n, rho, SNR, B, nbSimu, alpha)
  fTags <- c("B", "S", "Sbar", "V", "Vbar")
  filenames <- sprintf("Mein2006,%s,%s.xdr", fTags, simTags)
  pathnames <- file.path(path, filenames);
  mat <- matrix(NA, m, nbSimu)
  matList <- list()
  for (ff in seq(along=fTags)) {
    fTag <- fTags[ff]
    pathname <- pathnames[ff]
    if (file.exists(pathname)) {
      warning("File exists: ", pathname)
      if (!overwrite) {
        stop("Found existing results. Use 'overwrite=TRUE' to force re-calculation", )
      }
    }
    matList[[fTag]] <- mat
  }

  doit <- function(ii, verbose=FALSE) {
    if (verbose && (ii%%step==0)) {
      print(paste("Simulation", ii))
    }
    
    ## response
    y <- rbinom(n, 1, p)

    ## equi-correlated
    eps <- simulateGaussianNullsFromFactorModel(m, n=n, flavor="equi-correlated", rho=rho, cov=FALSE)

    ## means
    mu <- matrix(0, nrow=nrow(eps), ncol=ncol(eps)) ## m x n
    if (m0<m) {
      w1 <- which(y==1)
      mu[H1, w1] <- SNR*sqrt(100/n)
    }
    X <- mu+eps

    if (flavorH=="sansSouci") {
      ## NB:a single call to 'rowRanks' instead of B*m calls to 'rank'...
      w <- wilcoxStat(X, y, B=B)
      scoreMat <- w$stat0Mat
      o <- order(w$stat, decreasing=TRUE)
      stat <- w$stat[o]

      ## joint FWER control (through gammatification)
      if (stepDown) {
          resJ <- stepDownControl(stat, scoreMat, tau="kFWER", alpha=alpha, verbose=verbose)
      } else {
          resJ <- getJointFWERThresholds(scoreMat, tau="kFWER", alpha=alpha, maxSteps=1000)
      }
      thr <- resJ$thr
      BB <- sapply(stat, function(x) sum(x<=thr))  ## Eqn (7) in Meinshausen (2006) (*not* Vbar)
      ##      R <- sapply(thr, function(x) sum(stat>=x))   ## Number of rejections
      R <- 1:m
      Sbar <- pmax(0, cummax(R-BB[R]))
    } else {
      ## A faster test function for Wilcoxon's test
      testFUN <- function(x, y, ...) {
        mywilcox.test(x, y, p.value=FALSE, bypassTiesMethod=TRUE)
      }

      res <- howmany_dependent(t(X), y, test=testFUN, n.permutation=B, flavor=flavorH, sort=sort)
      o <- res$order
      BB <- res$boundingfunction   ## (loose) upper bound on number of FALSE discoveries among first rejections
      Sbar <- lowerbound(res)      ## lower bound on number of TRUE discoveries among first rejections

      ## sanity check
      R <- 1:m
      Sbar2 <- pmax(0, cummax(R-BB[R]))
      stopifnot(identical(Sbar, Sbar2))
      ## /sanity check
    }
    S <- cumsum(sapply(o, "%in%", H1))  ## number of TRUE discoveries among first rejections
    V <- cumsum(sapply(o, "%in%", H0))  ## number of FALSE discoveries among first rejections
    R <- S+V
    Vbar <- R-Sbar[R]                   ## (tighter) upper bound on number of FALSE discoveries among first rejections
    stopifnot(all(Vbar<=BB))  ## sanity check

    cbind(B=BB, S, Sbar, V, Vbar)
  }

  ## t0 <- Sys.time()
  ## for (kk in 1:nbSimu) {
  ##   res <- doit(kk)
  ##   matList[["BB"]][, kk] <- res[["BB"]]
  ##   matList[["S"]][, kk] <- res[["S"]]
  ##   matList[["Sbar"]][, kk] <- res[["Sbar"]] <- Sbar
  ##   matList[["V"]][, kk] <- res[["V"]]
  ##   matList[["Vbar"]][, kk] <- res[["Vbar"]]
  ## }

  resL <- parallel::mclapply(1:nbSimu, doit, verbose=verbose, mc.cores=mc.cores)
  resA <- simplify2array(resL)
  nms <- dimnames(resA)[[2]]
  stopifnot(identical(sort(fTags), sort(nms)))
  
  for (ff in seq(along=nms)) {
    fTag <- nms[ff]
    pathname <- pathnames[ff]
    saveObject(resA[, fTag, ], file=pathname)
  }
  res <- resA
}
