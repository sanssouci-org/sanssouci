<!-- badges: start -->
  [![R build status](https://github.com/pneuvial/sanssouci/workflows/R-CMD-check/badge.svg)](https://github.com/pneuvial/sanssouci/actions)
[![Coverage Status](https://codecov.io/gh/pneuvial/sanssouci/branch/develop/graph/badge.svg)](https://codecov.io/github/pneuvial/sanssouci/branch/develop)
 <!-- badges: end -->
 
The goal of sansSouci is to perform *post hoc inference*: in a multiple testing context, sansSouci provides statistical guarantees on possibly user-defined and/or data-driven sets of hypotheses. 


# Installation

```r
remotes::install_github("pneuvial/sanssouci@develop")
```

# Example use cases in genomics and neuroimaging

Typical use cases covered by sansSouci are:

- [differential gene expression studies](https://pneuvial.github.io/sanssouci/articles/post-hoc_differential-expression.html) in genomics
- [fMRI studies](https://pneuvial.github.io/sanssouci/articles/post-hoc_fMRI.html) in neuroimaging. 

In both cases, permutation-based post hoc inference typically outperforms classical post hoc bounds based on probabilistic inequalities.

sansSouci is developed within the [SansSouci project](https://www.math.univ-toulouse.fr/~pneuvial/sanssouci).