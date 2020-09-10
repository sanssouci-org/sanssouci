# kFWERThresholdFamily
#
# calibration of the 'balanced' threshold family of BNR
#
# This family \eqn{(t_k)} of thresholds is calibrated in such a way that for
# each \eqn{k\leq kMax}, \eqn{(t_k)} controls the (marginal) k-FWER
#
# The \eqn{m} x \eqn{B} matrix resulting of sorting the argument \code{mat}
# by row and then by columns is returned as the attribute 'Q' of the return
# value.
#
# @param mat A \eqn{m} x \eqn{B} matrix of Monte-Carlo samples of test
# statistics under the null hypothesis. \describe{ \item{m}{is the number of
# tested hypotheses} \item{B}{is the number of Monte-Carlo samples}}
# @param kMax For simultaneous control of (\eqn{k}-FWER for all \eqn{k \leq
# k[max]}).
# @param Rcpp If \code{TRUE} (the default), some costly operations (sorting)
# are performed in C++.
# @param verbose If \code{TRUE}, print results of intermediate calculations.
# Defaults to \code{FALSE}.
# @return A function \eqn{f:alpha \mapsto f(alpha)} such that for each
# \eqn{alpha}, \eqn{f(alpha)} is a candidate threshold family for controlling
# the joint FWER.
# @author Gilles Blanchard, Pierre Neuvial and Etienne Roquain
# @export
# @examples
#
# m <- 1023
# B <- 1e3
#
# flavor <- c("independent", "equi-correlated", "3-factor model")[2]
# rho <- 0.2
# mat <- gaussianTestStatistics(m, B, dep = "equi", param = rho)$X0
# pval <- pnorm(mat, lower.tail = FALSE)
# tk <- kFWERThresholdFamily(pval)
#
# thr <- tk(0.2)
#
# @export
#
kFWERThresholdFamily <- function(mat,
                                 kMax=nrow(mat),
                                 Rcpp=TRUE,
                                 verbose=FALSE) {
    Q <- -bisort(-mat, kMax=kMax, Rcpp=Rcpp)
    tk <- function(alpha) thresholdFamily(alpha, Q)
    attr(tk, 'Q') <- Q
    return(tk)
}
