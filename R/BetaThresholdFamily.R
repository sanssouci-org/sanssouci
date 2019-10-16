#' BetaThresholdFamily
#' 
#' Beta threshold family
#' 
#' @param m The number of hypotheses tested
#' @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \le 
#'   k[max]}).
#' @return A threshold function (on the scale of test statistics) which
#'   achieves "balanced" joint error rate (JER) control under independence
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @importFrom stats qbeta
BetaThresholdFamily <- function(m, kMax = m){
    sk <- function(alpha) qnorm(1-qbeta(alpha, 1:kMax, m+1-1:kMax))
    return(sk)
}
