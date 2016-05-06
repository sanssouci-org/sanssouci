# Package: sansSouci

## Version: 0.2.3 [2016-05-06]

* BUG FIX: in 'getJointFWERThresholds', function 'refFamily' could
  return a matrix of dimensions m x 0 instead of a vector of length m.

## Version: 0.2.2 [2016-04-20]

* Updated examples to pass R CMD check without ERROR.
* 'stepDownControl' vignette now in .Rmd.

## Version: 0.2.1 [2016-02-04]

* Fixed NAMESPACE for calling Rcpp.

## Version: 0.2.0 [2016-02-02]

* SPEEDUP: 'getJointFWERThresholds' now calls Rcpp functions for time costly operations (sorting).
* Fixed some typos in the documentation.

## Version: 0.1.3 [2016-01-28]

* BUG FIX: the one-parameter family should not be updated at each step down!
* Example fixed accordingly.
* Changed names of input and return values of 'getJointFWERThresholds' for consistency with the notation of the BNR paper.
* Added test script for unbalancedness.

## Version: 0.1.2 [2016-01-07]

* Implemented the "Oracle" version of step-down control.

## Version: 0.1.1 [2015-12-17]

* added scripts to test step-down JFWER control

## Version: 0.1.0 [2015-12-10]

* Tentative implementation of Etienne's trick to avoid dichotomy in 'getJointFWERThresholds'. Quite slow for kFWER thresholds.



