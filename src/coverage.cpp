#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector coverage(NumericVector thr, arma::mat Z) {
    // Calculate P_n(k-max(thr) > thr(k))
    int m = thr.size();
    int B = Z.n_cols;
    int m2 = Z.n_rows;
    if (m != m2) {
        stop("Length of argument 'thr' should match number of rows of argument 'Z'");
    }
    NumericVector prob;
    int count=0;
    for (int bb=0; bb<B; bb++) {
        for (int ii=0; ii<m; ii++) {
            if (Z(ii,bb)> thr(ii)) {
                count = count + 1;
                break;
            }
        }
    }
    prob = count/(double)B;
    return prob;
}
/*** R
m <- 10
B <- 20
thr <- sort(rnorm(m), decreasing=TRUE)
Z <- matrix(rnorm(m*B), m, B)
coverage(thr, Z)
*/
