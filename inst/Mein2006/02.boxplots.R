mat <- matrix(NA, nrow=length(rhos), ncol=length(ns))
colnames(mat) <- ns
rownames(mat) <- rhos
matM <- mat
matS <- mat

plot <- TRUE

for (nn in seq(along=ns)) {
  n <- ns[nn]
  print(paste("n=", n, sep=""))
  for (rr in seq(along=rhos)) {
    rho <- rhos[rr]
    print(paste("rho=", rho, sep=""))
    mvio <- boxplotsV(m, pi0, n, rho, nbSimu, path=path, figPath=figPath, plot=plot)
    if (!is.null(mvio)) {
      matM[rr, nn] <- mvio$mean
      matS[rr, nn] <- mvio$sd
    }
  }
}
