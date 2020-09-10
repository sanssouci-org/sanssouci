#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Sorting of the rows of a matrix by descending order
//
// @param X A numeric matrix.
// @return A numeric matrix whose rows are sorted by descending order.
//
// @export
// [[Rcpp::export]]
arma::mat rowSortDesc(arma::mat X) {
    // sorting by row in descending order
    int k = X.n_rows;
    int B = X.n_cols;
    arma::mat Y(k, B, arma::fill::zeros);
    for (int rr=0; rr<k; rr++) {
        arma::rowvec x = -X.row(rr);  // descending order
        std::sort(x.begin(), x.end());
        Y.row(rr) = -x;
    }
    return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
B <- rowSortDesc(A);
print(A);
*/
