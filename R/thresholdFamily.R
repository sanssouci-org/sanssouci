# thresholdFamily
#
# @param alpha A numeric value corresponding to \eqn{lambda(\alpha)} in BNR.
# @param Q A \eqn{m} x \eqn{B} matrix of \code{B} Monte-Carlo samples of 
# \code{m} p-values under the null hypothesis, sorted increasingly 
# by rows and then by columns.
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
# @export
thresholdFamily <- function(alpha, Q) {
    B <- ncol(Q)
    idx <- floor(min(alpha,1)*(B-1))
    if (idx == 0) {
        thr <- rep(0, nrow(Q))
    } else {
        thr <- Q[, idx]
    }
    thr
}
