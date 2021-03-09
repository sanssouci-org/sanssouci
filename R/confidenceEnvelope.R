#' Confidence bounds on the true/false positives among most significant items
#'
#' @param p.values A vector containing all \eqn{m} p-values, sorted
#'   increasingly
#'   
#' @param refFamily A character value, the reference family to be used. Should
#'   be either "Simes" (or equivalenlty, "Linear"), "Beta", or "Oracle".
#'   
#' @param param A numeric value or vector of parameters for the reference family. 
#' 
#' @param K For JER control over \code{1:K}, ie joint control of all
#'   \eqn{k}-FWER, \eqn{k \le K}. Defaults to \eqn{m}.
#' 
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#' 
#' \describe{
#' \item{FP}{Upper bound on the number of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the number of true positives in the 'x' most significant items}
#' \item{FDP}{Upper bound on the proportion of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the proportion of true positives in the 'x' most significant items}}
#' Defaults to \code{c("TP", "FDP")}.
#' 
#' @details \code{param} should be a numeric value unless \code{refFamily ==
#'   "Oracle"}. In the latter case, \code{param} should be a boolean vector of
#'   length \eqn{m} indicating whether each null hypothesis is true or false.
#'
#' @return A \code{data.frame} with \eqn{m} rows and 5 columns: \describe{
#' \item{x}{Number of most significant items selected}
#' \item{family}{Matches input argument \code{refFamily}}
#' \item{param}{Matches argument \code{param}}
#' \item{procedure}{Label for the procedure, typically of the form '<refFamily>(<param>)'}
#' \item{bound}{Value of the post hoc bound}
#' \item{stat}{Type of post hoc bound, as specified by argument \code{bound}}
#' }
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' rwt <- rowWelchTests(sim$X, sim$categ, alternative = "greater")
#' 
#' # calculate, print, and plot confidence bound
#' cb <- confCurveFromFam(rwt$p.value, refFamily = "Simes", param = 0.1)
#' head(cb)
#' plotConfCurve(cb, xmax = 200) 
#' 
confCurveFromFam <- function(p.values, refFamily, param, K = length(p.values), what = c("TP", "FDP")) {
    m <- length(p.values)
    fam0 <- c("Simes", "Beta", "Oracle")
    if (!(refFamily %in% fam0)) {
        stop("Unknown family: ", refFamily, "\n",
             "Only the following reference families are currently supported: ", 
             paste(fam0, collapse = ", "))
    }
    thr <- NULL
    if (refFamily %in% c("Simes", "Linear")) {
        thr <- SimesThresholdFamily(m, kMax = K)(param)
    } else if (refFamily == "Beta") {
        thr <- BetaThresholdFamily(m, kMax = K)(param)
    } else if (refFamily == "Oracle") {
        stopifnot(length(param) == m && all(param %in% c(0,1)))
        thr <- param
    }
    proc <- sprintf("%s(%s)", refFamily, param)
    bound(p.values, S = 1:m, thr, lab = proc, what = what, all = TRUE)
}



