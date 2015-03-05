## Retrieving results
flavor <- c("2006", "2013")[1]
study <- sprintf("Mein%s", flavor)

path <- file.path("resData", study)
path <- Arguments$getReadablePath(path)

filenames <- list.files(path)
pathnames <- file.path(path, filenames)

figPath <- file.path("fig", study, "fig1")
figPath <- Arguments$getWritablePath(figPath)

nbOver <- function(m, pi0, n, rho, nbSimu) {
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
  dep <- apply(Sbar>S, 2, sum)
  sum(dep>0)
}

nbSimu <- 1000
m <- 1000
pi0 <- 0.7
rhos <- c(0, 0.4)
ns <- c(20, 40, 60, 80)

tab <- matrix(NA, nrow=length(rhos), ncol=length(ns))
for (nn in seq(along=ns)) {
  for (rr in seq(along=rhos)) {
    tab[rr, nn] <- nbOver(m, pi0, ns[nn], rhos[rr], nbSimu)
  }
}
print(tab)
print(tab/nbSimu)

