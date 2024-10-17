#' Mini batch version of a rowTestFUN function.
#'
#' @param rowTestFUN a function taking as inputs Y, categ and Alternative and 
#' giving as a result a matrix of $p$-values. Example \code{\link{rowWelchTests}} or 
#' \code{\link{rowWilcoxonTests}}.
#' @param Y a $m \times n$ numeric matrix whose rows correspond to variables
#' and columns to observations
#' @param categ 
#' @param alternative 
#' @param max_batch_size 
#'
#' @return
#' @export
#'
#' @examples
mini_batch_rowTestFUN <- function(rowTestFUN, Y, categ,
                                  alternative = c("two.sided", "less", "greater"), 
                                  max_batch_size = 1e6){
  alternative <- match.arg(alternative)
  categ <- as.matrix(categ)

  m <- dim(Y)[1]
  B <- dim(categ)[2]
  nb_batch <- ceiling(B*m/max_batch_size)
  
  idxs <- rep(1:nb_batch, each = ceiling(B/nb_batch))[1:B]
  
  p_values <- matrix(NA_real_, nrow = m, ncol = B)
  for (batch_id in unique(idxs)) {
    id_batch <- which(idxs == batch_id)
    categ_batch <- categ[, id_batch]
    # fwt <- rowWelchtests.local(mat, categ_batch, alternative = "greater")
    p_values[, id_batch] <- rowTestFUN(Y, 
                                       categ_batch,
                                       alternative = alternative)$p.value
    # rm(fwt)
    # gc()
  }
  return(p_values)
}
