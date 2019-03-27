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
#' \item{bourgon}{Leukemia study of Chiaretti S, *et al.* (2004, 2005)}
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


#' Sample rMAPS data
#'
#' @format An integer matrix with 250 rows and 4642 columns. Each row
#'   corresponds to a genomic position relative to a splicing event, and each
#'   column corresponds to an exon. The entries of the matrix correspond to the
#'   number of matches of a particular DNA motif in the The samples are in two
#'   groups indicated by the column names of the matrix: 4550 controls ("0") and
#'   92 inclusion events ("1") corresponding to the status of the exon in a
#'   preliminary differential splicing analysis performed by rMATS.
#'   
#' @details rMAPS is a computational motif enrichment analysis tool for
#'   RNA-binding proteins using RNA-seq and CLIP-seq data. See
#'   http://rmaps.cecsresearch.org
#'   
#' @references Shen S., *et al.* (2014). rMATS: Robust and Flexible Detection of
#'   Differential Alternative Splicing from Replicate RNA-Seq Data. *PNAS*,
#'   111(51):E5593-601
#'   
#' @references Park JW, *et al.* (2016). rMAPS: RNA Map Analysis and Plotting
#'   Server for Alternative Exon Regulation. *Nucleic Acids Research*.
#' 
"rMAPS"

