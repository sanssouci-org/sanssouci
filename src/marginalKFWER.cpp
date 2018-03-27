#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Marginal Empirical coverage of a threshold family
//
// @param thr A numeric vector of length \code{m}, a threshold family. Should
//   be sorted in decreasing order.
// @param Z A \code{m} x \code{B} matrix of \code{B} Monte-Carlo samples of test statistics
//   under the null hypothesis. Each column should be sorted in decreasing order.
// @return The marginal empirical coverage of \code{thr} according to \code{Z}.
//
// @details The marginal empirical coverage of \code{thr} according to
//   \code{Z} is defined for an hypothesis of index \code{k} in \eqn{1 \dots m}
//   as the empirical probability the \eqn{k}-th max of the joint distribution of
//   the test statistics under the null hypothesis exceeds the \eqn{k}-th
//   largest value of \code{thr}, that is, \code{thr[k]}.
//
// @export
// [[Rcpp::export]]
NumericVector marginalKFWER(NumericVector thr, arma::mat Z) {
    // Calculate P_n(k-max(thr) > thr(k))
    int m = thr.size();
    int B = Z.n_cols;
    int m2 = Z.n_rows;
    if (m != m2) {
        stop("Length of argument 'thr' should match number of rows of argument 'Z'");
    }
    NumericVector probs(m);
    int count;
    for (int ii=0; ii<m; ii++) {
        count=0;
        double val = thr(ii);
        arma::rowvec z = Z.row(ii);
        for (int bb=0; bb<B; bb++) {
            if (z(bb)> val) {
             count = count+1;
            }
        }
        probs(ii) = count/(double)B;
    }
    return probs;
}

/*** R
m <- 10
B <- 20
thr <- sort(rnorm(m), decreasing=TRUE)
Z <- matrix(rnorm(m*B), m, B)
marginalKFWER(thr, Z)
*/
