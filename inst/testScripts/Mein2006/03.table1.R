if (FALSE) {  ## if run on its own: make sure everything is defined
  library(R.utils)
  sourceDirectory("package/inst/testScripts/Mein2006/R")
  study <- "Mein2006,sansSouci,0.0.4"
  path <- file.path("resData", study)
  path <- Arguments$getReadablePath(path)

  m <- 1e3
  nbSimu <- 1e4
  ns <- c(20, 60, 100)
  alpha <- 0.05

  rhos <- c(0, 0.2, 0.4)
  pi0 <- c(0.6, 0.99, 1)[3]
  SNR <- 1
}

mat <- matrix(NA, nrow=length(rhos), ncol=length(ns))
colnames(mat) <- ns
rownames(mat) <- rhos
matM <- mat
matS <- mat

for (nn in seq(along=ns)) {
  n <- ns[nn]
  print(paste("n=", n, sep=""))
  for (rr in seq(along=rhos)) {
    rho <- rhos[rr]
    print(paste("rho=", rho, sep=""))
    vio <- violators(m=m, pi0=pi0, n=n, rho=rho, B=B, nbSimu=nbSimu, alpha=alpha, SNR=SNR, path=path)
    if (!is.null(vio)) {
      matM[rr, nn] <- vio$mean
      matS[rr, nn] <- vio$sd
    }
  }
}

se <- round(100*2/sqrt(nbSimu)*matS, 1)
mat <- sprintf("$%s\\pm%s$", 100*matM, se)
dim(mat) <- dim(matM)
mat <- as.matrix(mat)
mat[is.na(matM)] <- NA_character_
colnames(mat) <- sprintf("$n=%s$", ns)
rownames(mat) <- sprintf("$\\rho=%s$", rhos)
mat

cap <- paste("The probability (in \\%) of underestimating the true proportion of false discoveries (using $\\alpha=", alpha, "$) for some $t\\in[0,1]$; $m_1=", round(m*(1-pi0)), "$.", sep="")
align <- paste(rep("l", ncol(mat)+1), collapse="")

library(xtable)
print(xtable(mat, caption=cap, align=align), type = "latex", sanitize.text.function = function(x){x})
