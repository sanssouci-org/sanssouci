#' post hoc bound obtained from Simes' inequality
#'
#' Lower bound on the number of correct rejections using Simes' reference
#' family
#'
#' If (R_k)_k provides jFWER control at level \eqn{\alpha} then with
#' probability greater than \eqn{1-\alpha}, \eqn{|H_0 cap R| \leq \min_k {|R
#' \cap (R_k)^c|+k-1}} A bit better: \eqn{|H_0 cap R| \leq (\min_{k<= |R|} {|R
#' \cap R_k^c|+k-1}) \wedge |R|}
#'
#' @param p A numeric vector of \code{m} p-values for all tested hypotheses.
#' @param select A vector of indices in \eqn{[1, \dots m]} of the hypotheses to be selected.
#' @param alpha A numeric value, the significance level of the test procedure.
#' @param Rcpp A boolean value: use Rcpp version of the implementation? Defaults to FALSE.
#' @param verbose If \code{TRUE}, prints verbose result to the
#'  screen. Defaults to \code{FALSE}.
#' @return A integer value, Simes's lower bound on the number of
#'  correct rejections within the selected hypotheses
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @references Goeman, J. J., & Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4), 584-597.
#' @export
#' @describeIn posthocBySimes R version
#' 
#' @examples
#' m <- 1e3
#' m1 <- 200
#' p <- 1-pnorm(c(rnorm(m1, mean=4), rnorm(m-m1, mean=0)))
#' R <- union(1:10, sample(m, 10))
#' alpha <- 0.10
#' if (require("cherry")) {
#'   hom <- hommelFast(p)
#'   pickSimes(hom, R, silent=TRUE, alpha = alpha)
#' }
#' posthocBySimes(p, R, alpha=alpha)
#'
posthocBySimes <- function(p, select, alpha, Rcpp = FALSE, verbose = FALSE) {
    m <- length(p)
    FUN <- posthocBySimes0
    if (Rcpp) {
        FUN <- posthocBySimes0Rcpp
    }
    m1hat <- FUN(p, 1:m, alpha = alpha)             # usual bound on total signals 
    pi0hat <- 1 - m1hat/m                             # corresponding bound on pi0
    FUN(p, select, alpha = alpha/pi0hat)
}

posthocBySimes0 <- function(p, select, alpha, verbose=FALSE) {
    m <- length(p)
    nR <- length(select)
    tSimes <- alpha*1:nR/m
    pR <- p[select]
    card <- sapply(tSimes, FUN=function(thr) {
                       sum(pR>thr)
                   })
    bounds <- pmin(card + 1:nR-1, nR)
    maxFalseRej <- min(bounds)
    minCorrRej <- nR - maxFalseRej
    minCorrRej
}
