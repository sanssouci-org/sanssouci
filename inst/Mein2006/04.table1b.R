library(R.utils)
sourceDirectory("package/inst/testScripts/Mein2006/R")
path <- "resData/Mein2006,sansSouci,0.0.5"

m <- 1e3
nbSimu <- 1e4
n <- 100
alpha <- 0.05

rhos <- c(0, 0.2, 0.4)
pi0s <- c(0.6, 0.99, 1)

mat <- matrix(NA, nrow=length(rhos), ncol=length(pi0s))
colnames(mat) <- pi0s
rownames(mat) <- rhos
matM <- mat
matS <- mat

for (pp in seq(along=pi0s)) {
  pi0 <- pi0s[pp]
  print(paste("pi0=", pi0, sep=""))
  for (rr in seq(along=rhos)) {
    rho <- rhos[rr]
    print(paste("rho=", rho, sep=""))
    vio <- violators(m, pi0, n, rho, B, nbSimu, alpha=alpha, path=path)
    if (!is.null(vio)) {
      matM[rr, pp] <- vio$mean/alpha
      matS[rr, pp] <- vio$sd/alpha
    }
  }
}

mea <- round(matM, 2)
se <- round(1/sqrt(nbSimu)*matS, 3)
mat <- sprintf("$%s\\pm%s$", mea, 2*se)
dim(mat) <- dim(matM)
mat <- as.matrix(mat)
mat[is.na(matM)] <- NA_character_
colnames(mat) <- sprintf("$\\pi_0=%s$", pi0s)
rownames(mat) <- sprintf("$\\rho=%s$", rhos)
mat

cap <- paste("$\\hat{p}/\\alpha$, where $\\alpha=", alpha, "$ and $\\hat{p}$ is the estimated probability of underestimating the true proportion of false discoveries for some $t\\in[0,1]$. $\\hat{p}$ is calculated based on $N=", nbSimu,  "$ simulation runs using the parameters $n=", n, "$  observations, $m=", m, "$ hypotheses, and $B=", B, "$ permutations.", sep="")
align <- paste(rep("l", ncol(mat)+1), collapse="")

library(xtable)
print(xtable(mat, caption=cap, align=align), type = "latex", sanitize.text.function = function(x){x})

tag <- sprintf("n=%s,alpha=%s", n, alpha)

## mean
what <- mea*alpha
filename <- sprintf("resM,%s.xdr", tag)
pathname <- file.path(path, filename)
saveObject(what, file=pathname)

## standard error
what <- se*alpha
filename <- sprintf("resS,%s.xdr", tag)
pathname <- file.path(path, filename)
saveObject(what, file=pathname)


