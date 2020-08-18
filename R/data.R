#' Sample data from differential gene expression studies
#'
#' A dataset containing microarray gene expression data from 4 cancer studies.
#'
#' @format A data frame with 21196 rows and 5 variables:
#' \describe{
#'   \item{id}{gene or probe set identified}
#'   \item{dataSet}{name of the study}
#'   \item{meanDiff}{difference between group means}
#'   \item{p.value}{two-sample Welch test p-value for the difference between group means}
#'   \item{adjp.BH}{FDR-adjusted p-values by the BH method}
#' }
#' 
#' @source The studies are the following: \describe{
#' \item{chiaretti}{Leukemia study of Chiaretti S, *et al.* (2004, 2005)}
#' \item{golub}{Leukemia study of Golub *et al.* (1999): 3051 genes and 38 tumor mRNA samples. Processed from \code{\link[multtest]{golub}}}
#' \item{hedenfalk}{Breast cancer study of Hedenfalk *et al.* (2001): 3226 genes and 5 samples. The two groups are BRCA1- and BRCA2-mutated patients.}
#' \item{rosenwald}{Diffuse large B-cell lymphoma of Alizadeh, A. A. *et al.* (2000). Processed from \code{\link[DLBCL]{exprLym}}}
#' }
#' 
#' 
#' @references Alizadeh, A. A. *et al.* (2000). Distinct types of diffuse large
#'   B-cell lymphoma identified by gene expression profiling. *Nature* 403:
#'   503-11.
#' 
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false
#'   discovery rate: a practical and powerful approach to multiple testing.
#'   *JRSS B*, 57(1), 289-300.
#' 
#' @references Chiaretti S, *et al.* (2004). Gene expression profile of adult
#'   T-cell acute lymphocytic leukemia identifies distinct subsets of patients
#'   with different response to therapy and survival. *Blood* 103:2771-2778.
#' 
#' @references Chiaretti S, *et al.* (2005). Gene expression profiles of
#'   B-lineage adult acute lymphocytic leukemia reveal genetic patterns that
#'   identify lineage derivation and distinct mechanisms of transformation.
#'   *Clinical Cancer Research* 11:7209-7219.
#' 
#' @references Golub, T. *et al.* (1999). Molecular classification of cancer:
#'   class discovery and class prediction by gene expression monitoring.
#'   *Science*, 286(5439), 531-537.
#'
#' @references Hedenfalk, I., *et al.* (2001). Gene-expression profiles in
#'   hereditary breast cancer. *New England Journal of Medicine*, 344(8),
#'   539-548.
#' 
"volcano"
