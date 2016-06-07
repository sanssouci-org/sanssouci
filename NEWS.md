# Package: sansSouci

## Version 0.4.1 [2016-06-07]

* BUG fix in 'stepDownJointFWERControl': the threshold family was not
  updated properly throughout the step-down process.

## Version 0.4.0 [2016-06-06]

Another major change in the package interface:
* added functions 'SimesPivotalStatistic', 'jointFWERThresholdCalibration', 'kFWERPivotalStatistic'.
* the main user-level entry function is now 'jointFWERControl', which performs step-down JFWER control by default.

## Version 0.3.2 (in progress)
* continuing the transition from inlienedocs to roxygen2 for package documentation.

## Version 0.3.1 [2016-06-02]
* BUG FIX in 'stepDownJointFWERControl': in some rare situations the algorithm could be stuck alternating between two candidate rejection sets R1, hence never terminating. Now forcing the 'best' candidate (ie the largest one) to be selected.
* Added a test (using 'extdata/simTrace.rda') to make sure that 'stepDownJointFWERControl' terminates in such situations.

## Version: 0.3.0 [2016-06-02]

Major rewriting of the code for increased clarity and efficiency:
* Renamed and reshaped the functions to calibrate JFWER thresholds (step down and single step).
* Dropped dichotomy in favor of pivotal statistic for the computation of JFWER thresholds.
* Speedup in 'upperBoundFP'.
* Added 'SimesThresholdFamily', 'pivotalStat', 'thresholdFamily', 'kFWERThresholdFamily', 'partialColSortDesc', 'biSort'.
* Added test scripts (using package "testthat").
* Moved from "inlinedocs" to "roxygen2" for documentation (mostly because documenting Rcpp functions is easier and the nice integration with other (dev)tools for checking package, running examples and performing package tests).

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



