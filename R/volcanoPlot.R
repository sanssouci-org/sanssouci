#' Volcano plot
#' 
#' Volcano plot for differential expression studies
#' 
#' @param dat A numeric matrix whose rows correspond to variables and columns to
#'   observations

#' @param thr A numeric vector of length K, a JER controlling family
#' 
#' @param categ A vector of \code{ncol(mat)} categories in \eqn{'0','1'} for the
#'   observations

#' @param p A numeric value, the p-value threshold under which genes are selected
#' @param q A numeric value, the q-value (or FDR-adjusted p-value) threshold under which genes are selected
#' @param r A numeric value, the absolute fold change above which genes are selected
#' @param cex A numeric vector of length 2, the relative magnification factor for unselected (\code{cex[1]}) and unselected (\code{cex[2]}) genes. 
#' 
#' @param col A vector of length 3, 
#' @param pch An integer or single character string specifying the plotting character, see \code{\link{par}}
#' @param ylim A numeric vector of length 2, the \eqn{y} limits of the plot
#'
#' @details A Welch T-test of differential expression between the two categories
#'   defined by \code{categ} are applied for each gene using the
#'   \code{\link{rowWelchTests}} function, which also outputs the "fold change"
#'   (mean difference in log scale) between the two categories.
#'   
#' @export
#' @importFrom graphics abline legend rect title
#' @importFrom stats p.adjust
#' 
#' @examples
#' m <- 500
#' pi0 <- 0.5
#' m1 <- m-m*pi0
#' sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
#'                        pi0 = pi0, SNR = runif(m1)*6-3, prob = 0.5)
#' X <- sim$X
#' alpha <- 0.2
#' cal <- calibrateJER(X, B = 1e2, alpha = alpha, refFamily="Simes")
#' volcanoPlot(X, categ = colnames(X), thr = cal$thr, p = 3, r = r, ylim = c(0, 6))


volcanoPlot <- function(dat, thr, categ = colnames(dat), 
                        p = 1, q = 1, r = 0,
                        cex = c(0.2, 0.6), 
                        col = c("#33333333", "#FF0000", "#FF666633"),
                        pch = 19, ylim = NULL) {
    if (p <1 && q < 1) {
        warning("Filtering both on p-values and BH-adjusted p-values")
    }
    m <- nrow(dat)
    
    ## Student/Welch tests
    dex <- rowWelchTests(dat, categ)

    ## p-values
    pval <- dex[["p.value"]]
    logp <- -log10(pval)

    ## fold changes
    fc <- dex$meanDiff  
    
    adjp <- p.adjust(pval, method = "BH")  ## adjusted p-values
    bhsel <- (adjp <= q)                   ## selected by BH at level q
    bh <- min(logp[bhsel])                 ## threshold on the p-value scale
    
    ## gene selections
    sel1 <- which(logp >= bh & fc >= r)
    sel2 <- which(logp >= bh & fc <= -r)
    sel12 <- union(sel1, sel2)
    
    ## post hoc bounds in selections
    n1 <- length(sel1)
    FP1 <- maxFP(pval[sel1], thr = thr)
    TP1 <- n1 - FP1
    FDP1 <- round(FP1/n1, 2)
    
    n2 <- length(sel2)
    FP2 <- maxFP(pval[sel2], thr = thr)
    TP2 <- n2 - FP2
    FDP2 <- round(FP2/n2, 2)
    
    n12 <- length(sel12)
    FP12 <- maxFP(pval[sel12], thr = thr)
    TP12 <- n12 - FP12
    FDP12 <- round(FP12/n12, 2)
    
    ## graphical parameters
    cols <- rep(col[1], m)
    cols[c(sel1, sel2)] <- col[2]
    
    cexs <- rep(cex[1], m)
    cexs[sel12] <- cex[2]

    xlab <- "Fold change (log scale)"
    ylab <- expression("p-value (-log[10] scale)")
    infty <- 100
    if (is.null(ylim)) {
        ylim <- c(0, max(logp))
    }
    plot(fc, logp, pch = pch, cex = cexs, col = cols, 
         xlab = xlab, ylab = ylab, ylim = ylim)
    rect(xleft = -infty, ybottom = bh, xright = -r, ytop = infty, 
         col = col[3], border = NA, lwd = 2)
    rect(xleft = r, ybottom = bh, xright = infty, ytop = infty, 
         col = col[3], border = NA, lwd = 2)
    abline(h = bh, col = "gray")
    abline(v = c(-1, 1)*r, col = "gray")
    txt <- c(sprintf("%s genes\nTP > %s\nFDP < %s", n2, TP2, FDP2))
    legend("topleft", txt, border = "white", bty = "n", text.col = col[2])
    
    txt <- c(sprintf("%s genes\nTP > %s\nFDP < %s", n1, TP1, FDP1))
    legend("topright", txt, border = "white", bty = "n", text.col = col[2])
    
    txt <- c(sprintf("%s genes selected\nAt least %s true positives (FDP < %s)", 
                     n12, TP12, FDP12))
    title(txt)
}