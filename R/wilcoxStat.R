##' wilcoxStat
##'
##' Calculation of Wilcoxon sum rank test statistic
##'
##'
##' @param X An \eqn{m x n} covariate matrix
##' @param y A vector of \eqn{n} phenotypes in \eqn{{0,1}}
##' @param B Number of resamplings
##' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
##' @return A list with elements \describe{
##'   \item{stat}{A vector of \code{m} Wilcoxon sum rank test statistics of association between \code{X} and \code{y}.}
##'   \item{stat0Mat}{An \code{m} x \code{B} matrix of \code{B} realizations of a \code{m}-dimensional vector of test statistics under the null hypothesis of no association between \code{X} and \code{y}.}}
##' @importFrom matrixStats rowRanks
##' @export
wilcoxStat <- function(
    ### Fast calculation of Wilcoxon sum rank test statistic
    X,
    ### An \eqn{m x n} covariate matrix
    y,
    ### A vector of \eqn{n} phenotypes in \eqn{{0,1}}
    B=1000
    ### Number of resamplings
    ) {
  rks <- rowRanks(X)
  n0 <- sum(y==0)
  n1 <- sum(y==1)
  rs <- n0*(n0+1)/2  ## mininmal rank sum
  med <- n0*n1/2

  w0 <- which(y==0)
  stat <- rowSums(rks[, w0])

  stat0 <- replicate(B, {
    yy <- sample(y)
    w0 <- which(yy==0)
    rowSums(rks[, w0])
  })

  list(
      stat=abs(stat-rs-med),
      ### A vector of \eqn{m}: Wilcoxon sum rank test statistics of association between \code{X} and \code{y}
      stat0Mat=abs(stat0-rs-med)
      ### A \eqn{m x B} matrix: \eqn{B} realizations of a \eqn{m}-dimensional vector of test statistics under the null hypothesis of no association between \code{X} and \code{y}.
      )
}
