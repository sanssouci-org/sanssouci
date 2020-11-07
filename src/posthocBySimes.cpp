#include <Rcpp.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <cmath>
#include <vector>
using namespace Rcpp;

NumericVector pmin2(NumericVector x, double y) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        out[i] = std::min(x[i], y);
    }
    
    return out;
}

//' @describeIn posthocBySimes Rcpp version
//' @export
// [[Rcpp::export]]
double posthocBySimes0Rcpp(NumericVector p, NumericVector select, double alpha){
    double m = p.length();
    double nR = select.length();
    
    NumericVector Seq;
    for(int i=1; i<=nR ; i++)
    {
        Seq.push_back(i);
    }
    
    NumericVector tSimes(nR);
    for(int i = 0; i < nR; i++) tSimes[i] = Seq[i]*alpha/m;
    // transform(Seq.begin(), Seq.end(), Seq.begin(),
    //            bind2nd(std::multiplies<double>(),alpha/m));
    
    NumericVector pR = p[select-1];
    
    NumericVector card(tSimes.length());
    for(int i=0; i<tSimes.length(); i++){
        card[i] = sum(pR > tSimes[i]) + (Seq[i]-1);
    }
    NumericVector bounds = pmin2(card,nR);
    double maxFalseRej = *std::min_element(bounds.begin(), bounds.end());
    double minCorrRej = nR - maxFalseRej;
    
    return minCorrRej ;
}

/*** R
data(NAEP, package="cherry")
posthocBySimes0Rcpp(NAEP, 1:30, 0.05)
posthocBySimes0(NAEP, 1:30, 0.05)
*/

