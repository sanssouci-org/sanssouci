#' Confidence bounds on the true/false positives among most significant items
#'
#' @param p_values A vector containing the p-values for all m hypotheses,
#'   sorted increasingly
#'
#' @param ref_family A character value, the reference family to be used. Should
#'   be either "Simes" (or equivalenlty, "Linear"), "Beta", or "Oracle".
#'
#' @param param A numeric value or vector of parameters for the reference
#' family.
#'
#' @param K For JER control over `1:K`, ie joint control of all
#'   k-FWER, k <= K. Defaults to m.
#'
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#'
#' - FP: Upper bound on the number of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the number of true positives in the 'x' most
#' significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most
#' significant items.
#'
#' Defaults to `c("TP", "FDP")`
#'
#' @details `param` should be a numeric value unless `ref_family ==
#'   "Oracle"`. In the latter case, `param`` should be a boolean vector of
#'   length m indicating whether each null hypothesis is true or false.
#'
#' @return A `data.frame` with `m` rows and 5 columns:
#' - x: Number of most significant items selected
#' - family: Matches input argument `ref_family`
#' - param: Matches argument `param`
#' - procedure: Label for the procedure, typically of the form
#' 'ref_family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`
#'

#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @export
#' @examples
#'
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
#' prob = 0.5)
#' rwt <- row_welch_tests(sim$X, sim$categ, alternative = "greater")
#'
#' # calculate, print, and plot confidence bound
#' cb <- conf_curve_from_fam(rwt$p.value, ref_family = "Simes", param = 0.1)
#' head(cb)
#' plotConfCurve(cb, xmax = 200)
#'
conf_curve_from_fam <- function(p_values, ref_family, param,
                                K = length(p_values), what = c("TP", "FDP")) {
  m <- length(p_values)
  fam0 <- c("Simes", "Linear", "Beta", "Oracle")
  if (!(ref_family %in% fam0)) {
    stop(
      "Unknown family: ", ref_family, "\n",
      "Only the following reference families are currently supported: ",
      paste(fam0, collapse = ", ")
    )
  }
  thr <- NULL
  if (ref_family %in% c("Simes", "Linear")) {
    thr <- t_linear(param, seq_len(K), m)
  } else if (ref_family == "Beta") {
    thr <- t_linear(param, seq_len(K), m)
  } else if (ref_family == "Oracle") {
    stopifnot(length(param) == m && all(param %in% c(0, 1)))
    thr <- param
  }
  proc <- sprintf("%s(%s)", ref_family, param)
  posthoc_bound(p_values, S = 1:m, thr, lab = proc, what = what, all = TRUE)
}
