posthoc <- structure(function(
    ### Lower bound on the number of correct rejections.
    p, ### A vector of p-values for all tested hypotheses.
    select,  ### The indexing vector of the p-values of the hypotheses to be selected.
    alpha=0.05, ### The significance level of the test procedure.
    silent=FALSE ### If \code{FALSE}, prints verbose result to the screen.
    ) {
  ##details<<If (R_k)_k provides jFWER control at level \eqn{\alpha} then
  ##  \eqn{|H_0 cap R| \leq \min_k {|R \cap (R_k)^c|+k-1}}
  ## A bit better:
  ##  \eqn{|H_0 cap R| \leq (\min_{k<= |R|} {|R \cap R_k^c|+k-1}) \wedge |R|}
  m <- length(p)
  o <- order(p)
  po <- p[o]
  nR <- length(select)
  ## tSimes <- alpha*1:m/m
  tSimes <- alpha*1:nR/m
  pR <- p[select]
  card <- sapply(tSimes, FUN=function(thr) {
    sum(pR>thr)
  })
  ## bounds <- pmin(card + 1:m-1, m)
  bounds <- pmin(card + 1:nR-1, nR)
  maxFalseRej <- min(bounds)
  minCorrRej <- nR - maxFalseRej

  if (!silent) {
    cat("Rejected ", nR, " hypotheses. At confidence level ", 1 - alpha, ":\n", sep = "")
    cat("Correct rejections >= ", minCorrRej, "; ", sep = "")
    cat("False rejections <= ", maxFalseRej, ".\n", sep = "")
    invisible(minCorrRej)
  }
  else {
    minCorrRej
  }
### Lower bound on the number of correct rejections.
}, ex=function() {
  if (require(cherry)) {
    data("NAEP", package="cherry")
    p <- NAEP
    R <- c("HI","MN","IA")
  } else {
    m <- 1e2
    m0 <- 10
    p <- 1-pnorm(c(rnorm(m0, mean=3), rnorm(m-m0, mean=0)))
    R <- 1:15
  }
  posthoc(p, R)

  if (require(cherry)) {   ## Comparison with 'cherry picking'
    pickSimes(p, R)
  }
  
  if (require(cherry)) {   ## Comparison with 'cherry picking'
    kk <- 0
    res <- NULL
    diffIdxs <- NULL
    while (TRUE) {
      kk <- kk+1
      if (kk %% 1000 == 0) {
        cat(paste("Iteration: ", kk, "\n", sep=""))
      }
      m <- 1e3
      m0 <- 100
      p <- 1-pnorm(c(rnorm(m0, mean=4), rnorm(m-m0, mean=0)))
      nR <- min(m, rpois(1, lambda=100))
      R <- sample(1:m, nR)
      ph <- posthoc(p, R, silent=TRUE)
      ps <- pickSimes(p, R, silent=TRUE)
      res <- c(res, ph)
      diff <- ph-ps
      if (diff!=0) {
        diffIdxs <- c(diffIdxs, kk)
        print(diff)
        ## break;
      } 
      stopifnot(abs(diff)<=1)
    }
  }
})
