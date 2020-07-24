#' SimesThresholdFamily
#'
#' Simes' threshold family
#'
#' @param m The number of hypotheses tested
#' @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le
#' k[max]}).
#' @return A threshold function (on the scale of test statistics) based on the
#' classical family of thresholds introduced by Simes (1986):
#' \eqn{\alpha*k/m}. This family yields joint FWER control at level (at most)
#' \eqn{\alpha} if the test statistics are positively dependent (PRDS) under
#' H0.
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#'
#' sk <- SimesThresholdFamily(12)
#' thr <- sk(0.2)
#'
SimesThresholdFamily <- function(m, kMax = m){
    tk <- function(alpha) min(alpha, 1)*(1:kMax)/m
    return(tk)
}
