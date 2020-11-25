#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Sorting of the colums of a matrix by ascending order
//
// @param X A numeric matrix.
// @return A numeric matrix whose rows are sorted by ascending order.
//
// @export
// [[Rcpp::export]]
arma::mat colSort(arma::mat X) {
    int m = X.n_rows;
    int B = X.n_cols;
    arma::mat Y(m, B, arma::fill::zeros);
    for (int cc=0; cc<B; cc++) {
        arma::colvec x = X.col(cc); 
        std::sort(x.begin(), x.end());
        Y.col(cc) = x;
    }
    return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
B <- colSort(A);
print(A);
*/
