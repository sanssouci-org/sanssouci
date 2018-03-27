# Version 0.6.0 [2018-03-21]

## Several important user-level simplifications

* Added function maxFP and minTP to calculate post hoc bound among user-defined selections. Renamed function 'upperBoundFP' to 'curveMaxFP' for internal consistency. Removed functions 'falseDiscoveryEnvelope', 'falseDiscoveryProportionControl' and 'falseDiscoveryProportionEnvelope' which provided correct but conservative bounds. 
* Added function calibrateJER as the main user-level entry point for calibrating JER thresholds. 
* Cleanup: user-level entry point for simulating test statistics is now 'gaussianTestStatistics'.
* Cleanup: user-level entry point for simulating observations from one and two-sample problems is now 'gaussianSamples'.
* Added corresponding tests and updated vignettes accordingly.
* Cleanups in 'inst' scripts.
* Cleanup: *many* lower-level functions are not exported anymore; their documentation is still available in the form of comments but not as .Rd files. 
* Some code cleanups (spacing; comments)

# Version 0.5.1 [2018-02-07]

* Added data/volcano.rda with p-values and fold changes for 4 public gene expression data sets.
* Added first draft of a volcano plot shiny app.

# Version 0.5.0 [2017-12-20]

* Added/updated randomization tests (by permutation or sign-flipping).
* Completed numerical experiments for BNR paper.

# Version 0.4.12 [2017-10-28]

* Added stratified Monte-Carlo estimation of JER coverage in BNR experiments
  (as nicely suggested by a referee).
* Added simulateEquiByGroupPermutation.
* Added C++ verson of postHocBySimes (written by @sbsad) + corresponding test script.
* Added "twoGroupPermutationTest" (written by @sbsad) + corresponding test script.
* Added "rowWelchTest" + corresponding test script.

# Version 0.4.11 [2017-02-14]

* Updated step down vignette and documentation.
* Added travis CI and appveyor support.

# Version 0.4.10 [2017-01-27]

* Added tiny test for numerical reproducibility of JR control function.
* Speedup in 'kFWERPivotalStatistic': avoid calling 'partialColSortDesc' twice when 'C=H' (single step control)

# Version 0.4.9 [2017-01-10]

* Moved to github
* Passes R CMD check locally

# Version 0.4.8 [2016-08-11]

* Reshaped the scripts of the numerical experiments or the BNR paper.

# Version 0.4.7 [2016-07-04]

* In 'upperBoundFP', implemented the case 'kMax<m' for all flavors. Added corresponding 'testthat' scripts.

# Version 0.4.6 [2016-06-17]

* Added linearly decreasing SNR to BNR test scripts.
* Added 'simulateToeplitz' to simulate under Gaussian distribution
  with Toeplitz covariance matrix.

# Version 0.4.5 [2016-06-16]

* Added alternative (finer?) definition of power to BNR test scripts.

# Version 0.4.4 [2016-06-10]

* BUG fix in 'pivotatStat' (cont'd): the pivotal statistic was not
  calculated properly for kFWER threshold family when the cardinal of
  the set of candidate hypotheses ($\cC$) was less than m.
* Added a test script to check that the R and C++ versions of
  'minPseudoRanks' give identical results.

# Version 0.4.3 [2016-06-09]

* BUG fix in 'pivotatStat' (cont'd): the pivotal statistic was not
  calculated properly for kFWER threshold family when the cardinal of
  the set of candidate hypotheses ($\cC$) was less than m.
* Removed function 'pivotalStat' as the calculation of kFWER and Simes
  statistics now have almost nothing to share.

# Version 0.4.2 [2016-06-09]

* BUG fix in 'pivotatStat': the pivotal statistic was not calculated
  properly for Simes' threshold family when the cardinal of the set of
  candidate hypotheses ($\cC$) was less than m.
* removed 'jointFWERThresholdCalibration' to make the code in 'jointFWERControl' more explicit.

# Version 0.4.1 [2016-06-07]

* BUG fix in 'jointFWERControl': the threshold family was not
  updated properly throughout the step-down process.

# Version 0.4.0 [2016-06-06]

Another major change in the package interface:
* added functions 'SimesPivotalStatistic', 'jointFWERThresholdCalibration', 'kFWERPivotalStatistic'.
* the main user-level entry function is now 'jointFWERControl', which performs step-down JFWER control by default.

# Version 0.3.2 (in progress)
* continuing the transition from inlienedocs to roxygen2 for package documentation.

# Version 0.3.1 [2016-06-02]
* BUG FIX in 'stepDownJointFWERControl': in some rare situations the algorithm could be stuck alternating between two candidate rejection sets R1, hence never terminating. Now forcing the 'best' candidate (ie the largest one) to be selected.
* Added a test (using 'extdata/simTrace.rda') to make sure that 'stepDownJointFWERControl' terminates in such situations.

# Version: 0.3.0 [2016-06-02]

Major rewriting of the code for increased clarity and efficiency:
* Renamed and reshaped the functions to calibrate JFWER thresholds (step down and single step).
* Dropped dichotomy in favor of pivotal statistic for the computation of JFWER thresholds.
* Speedup in 'upperBoundFP'.
* Added 'SimesThresholdFamily', 'pivotalStat', 'thresholdFamily', 'kFWERThresholdFamily', 'partialColSortDesc', 'biSort'.
* Added test scripts (using package "testthat").
* Moved from "inlinedocs" to "roxygen2" for documentation (mostly because documenting Rcpp functions is easier and the nice integration with other (dev)tools for checking package, running examples and performing package tests).

# Version: 0.2.3 [2016-05-06]

* BUG FIX: in 'getJointFWERThresholds', function 'refFamily' could
  return a matrix of dimensions m x 0 instead of a vector of length m.

# Version: 0.2.2 [2016-04-20]

* Updated examples to pass R CMD check without ERROR.
* 'stepDownControl' vignette now in .Rmd.

# Version: 0.2.1 [2016-02-04]

* Fixed NAMESPACE for calling Rcpp.

# Version: 0.2.0 [2016-02-02]

* SPEEDUP: 'getJointFWERThresholds' now calls Rcpp functions for time costly operations (sorting).
* Fixed some typos in the documentation.

# Version: 0.1.3 [2016-01-28]

* BUG FIX: the one-parameter family should not be updated at each step down!
* Example fixed accordingly.
* Changed names of input and return values of 'getJointFWERThresholds' for consistency with the notation of the BNR paper.
* Added test script for unbalancedness.

# Version: 0.1.2 [2016-01-07]

* Implemented the "Oracle" version of step-down control.

# Version: 0.1.1 [2015-12-17]

* added scripts to test step-down JFWER control

# Version: 0.1.0 [2015-12-10]

* Tentative implementation of Etienne's trick to avoid dichotomy in 'getJointFWERThresholds'. Quite slow for kFWER thresholds.



