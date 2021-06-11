#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// partial sorting of the columns of a matrix by descending order
//
// @param X A numeric matrix
// @param k An integer value between 1 and \code{nrow(mat)}
// @export
// [[Rcpp::export]]
arma::mat partialColSortDescCpp(arma::mat X, int k) {
// partial sorting by column in descending order
    int B = X.n_cols;
    arma::mat Y(k, B, arma::fill::zeros);
    for (int bb=0; bb<B; bb++) {
        arma::vec x = -X.col(bb);  // descending order
        // std::partial_sort(x.begin(), x.begin()+k, x.end()); // very slow!
        std::nth_element(x.begin(), x.begin()+k, x.end());
        std::sort(x.begin(), x.begin()+k);
        for (int rr=0; rr<k; rr++) {
            Y(rr,bb) = -x(rr);    // back to original values
        }
    }
    return Y;
}

// [[Rcpp::depends(RcppArmadillo)]]

// partial sorting of the columns of a matrix by descending order
//
// @param X A numeric matrix
// @param k An integer value between 1 and \code{nrow(mat)}
// @export
// [[Rcpp::export]]
arma::mat partialColSortCpp(arma::mat X, int k) {
    // partial sorting by column
    int B = X.n_cols;
    arma::mat Y(k, B, arma::fill::zeros);
    for (int bb=0; bb<B; bb++) {
        arma::colvec x = X.col(bb);
        // std::partial_sort(x.begin(), x.begin()+k, x.end()); // very slow!
        std::nth_element(x.begin(), x.begin()+k, x.end());
        std::sort(x.begin(), x.begin()+k);
        Y.col(bb) = x.subvec(0, k-1);
    }
    return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
B <- partialColSortDescCpp(A, 2);
print(B);
C <- partialColSortCpp(A, 2);
print(C);

if (require("rbenchmark")) {
    A <- matrix(rnorm(1e6), 10000, 1000);
    kmax <- 1000;
    benchmark(partialColSortDescCpp(A, kmax),
              -apply(-A, 2, sort, partial=kmax), replications=1)
}
*/
