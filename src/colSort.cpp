#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Sorting of the colums of a matrix by ascending order
//
// @param X A numeric matrix.
// @return A numeric matrix whose columns are sorted by ascending order.
//
// @export
// [[Rcpp::export]]
arma::mat colSort(arma::mat X) {
    arma::mat Y = sort(X);
    return Y;
}

/*** R
A <- matrix(rnorm(15), 5, 3);
print(A)
apply(A, 2, sort)
colSort(A);

# A <- matrix(rnorm(1e7), 1e4, 1e3);
# b <- bench::mark(apply(A, 2, sort),
#                  colSort(A), iterations = 10)
# plot(b)
*/

