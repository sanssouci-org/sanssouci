#' post hoc bound obtained from Simes' inequality
#'
#' Lower bound on the number of correct rejections using Simes' reference
#' family
#'
#' @param p A numeric vector of \code{m} p-values for all tested hypotheses.
#' @param select A vector of indices in \eqn{[1, \dots m]} of the hypotheses to be selected.
#' @param alpha A numeric value, the significance level of the test procedure.
#' @return A integer value, Simes's lower bound on the number of
#'  correct rejections within the selected hypotheses
#' @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @references Goeman, J. J., & Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4), 584-597.
#' @export
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
#' p <- c(0, 0.01, 0.02, 0.05, 0.1)
#' alpha <- 0.2
#' R <- 1:length(p)
#' posthocBySimes(p, R, alpha, flavor = "single step")
#' posthocBySimes(p, R, alpha, flavor = "one step down")
#' posthocBySimes(p, R, alpha, flavor = "full step down")
#' if (require("cherry")) {
#'   hom <- hommelFast(p)
#'   pickSimes(hom, R, silent=TRUE, alpha = alpha)
#' }

posthocBySimes <- function(p, select, alpha, 
                           flavor = c("full step down", 
                                      "one step down",
                                      "single step")) {
  m <- length(p)
  stopifnot(all(select <= m))
  stopifnot(all(select >= 1))
  stopifnot(alpha >= 0)
  stopifnot(alpha <= 1)
  flavor <- match.arg(flavor)
  
  thr <- t_linear(alpha, 1:m, m)
  if (flavor == "one step down") {
    m0_hat <- maxFP(p.values = p, thr = thr) 
    thr <- t_linear(alpha, 1:m0_hat, m0_hat)
  } else if (flavor == "full step down") {
    m0_hat0 <- m
    m0_hat <- maxFP(p.values = p, thr = thr) 
    thr <- t_linear(alpha, 1:m0_hat, m0_hat)
    while (m0_hat < m0_hat0) {
      m0_hat0 <- m0_hat
      m0_hat <- maxFP(p.values = p, thr = thr) 
      thr <- t_linear(alpha, 1:m0_hat, m0_hat)
    }
  }
  minTP(p[select], thr = thr)
}

