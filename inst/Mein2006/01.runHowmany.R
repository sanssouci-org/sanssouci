for (n in ns) {
  for (rho in rhos) {
    st <- system.time(runMeinshausen2006(m, pi0, n, rho, B=B, nbSimu, alpha=alpha, SNR=SNR,
                                         flavorH=flavorH, path=path, overwrite=runForce, verbose=verbose, sort=sort, mc.cores=mc.cores))
    print(st)
    ## boxplotsV(m, pi0, n, rho, nbSimu, path=path, figPath=figPath)
    ## boxplots(m, pi0, n, rho, nbSimu, path=path, figPath=figPath)
  }
}
