#' Volcano plot
#' 
#' @param x An object. See individual methods for specifics
#' @param ... Other arguments passed to methods
#' @export
volcanoPlot <- function(x, ...) UseMethod("volcanoPlot")

#' @rdname volcanoPlot
#' @param x An object of class `SansSouci`
#' @param p A numeric value, the p-value threshold under which genes are selected
#' @param q A numeric value, the q-value (or FDR-adjusted p-value) threshold under which genes are selected
#' @param r A numeric value, the absolute fold change above which genes are selected
#' @param cex A numeric vector of length 2, the relative magnification factor for unselected (`cex[1]`) and unselected (`cex[2]`) genes. 
#' 
#' @param col A vector of length 3
#' @param pch An integer or single character string specifying the plotting character, see [par]
#' @param ylim A numeric vector of length 2, the `y` limits of the plot
#' @param ... Other arguments to be passed to volcanoPlot.numeric
#' @export
#' 
#' @examples
#' data(expr_ALL, package = "sansSouci.data")
#' groups <- ifelse(colnames(expr_ALL)=="NEG", 0, 1)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#' 
#' res <- fit(a, B = 100, alpha = 0.1)
#' volcanoPlot(res, q = 0.2, r = 0.2, ylim = c(0, 4))
volcanoPlot.SansSouci <- function(x, 
                                  p = 1, q = 1, r = 0,
                                  cex = c(0.2, 0.6), 
                                  col = c("#33333333", "#FF0000", "#FF666633"),
                                  pch = 19, ylim = NULL, ...) {
    object <- x; rm(x);
    if (object$input$n_group == 1) {
        stop("Can't do a volcano plot for one-sample tests!")
    }
    pval <- pValues(object)
    fc <- foldChanges(object)
    thr <- thresholds(object)
    volcanoPlot(pval, fc = fc, thr = thr, 
                p = p, q = q, r = r,
                cex = cex, 
                col = col,
                pch = pch, ylim = ylim, ...)
}

#' Volcano plot
#' 
#' Volcano plot for differential expression studies
#' 
#' @param x A vector of p-values
#' @param fc A vector of fold changes, of the same length as `x`
#' @param thr A numeric vector of length K, a JER controlling family
#' @param p A numeric value, the p-value threshold under which genes are selected
#' @param q A numeric value, the q-value (or FDR-adjusted p-value) threshold under which genes are selected
#' @param r A numeric value, the absolute fold change above which genes are selected
#' @param cex A numeric vector of length 2, the relative magnification factor for unselected (\code{cex[1]}) and unselected (\code{cex[2]}) genes. 
#' 
#' @param col A vector of length 3
#' @param pch An integer or single character string specifying the plotting character, see \code{\link{par}}
#' @param ylim A numeric vector of length 2, the \eqn{y} limits of the plot
#' @param bounds A boolean value: should the post hoc bounds be displayed on the plot? Defaults to TRUE
#' @param ... Not used
#'
#' @details A Welch T-test of differential expression between the two categories
#'   defined by \code{categ} are applied for each gene using the
#'   \code{\link{rowWelchTests}} function, which also outputs the "fold change"
#'   (mean difference in log scale) between the two categories.
#' @return The indices of selected genes (returned invisibly)
#' 
#' @importFrom graphics abline legend rect title
#' @importFrom stats p.adjust
#' @seealso Volcano plot shiny app at \url{ https://shiny-iidea-sanssouci.apps.math.cnrs.fr/}
volcanoPlot.numeric <- function(x, fc, thr,
                        p = 1, q = 1, r = 0,
                        cex = c(0.2, 0.6), 
                        col = c("#33333333", "#FF0000", "#FF666633"),
                        pch = 19, ylim = NULL, bounds = TRUE, ...) {
    pval <- x; rm(x);
    if (p < 1 && q < 1) {
        warning("Filtering both on p-values and BH-adjusted p-values")
    }
    m <- length(pval)
    
    ## sanity checks
    stopifnot(length(fc) == m)
    stopifnot(length(thr) <= m)
    
    logp <- -log10(pval)
    adjp <- p.adjust(pval, method = "BH")  ## adjusted p-values
    y_sel <- which((adjp <= q) &           ## selected by q-value
                       (pval <= p))        ##          or p-value
    y_thr <- Inf
    if (length(y_sel) > 0) {
        y_thr <- min(logp[y_sel])       ## threshold on the log(p-value) scale
    }
    
    ## gene selections
    sel1 <- which(logp >= y_thr & fc >= r)
    sel2 <- which(logp >= y_thr & fc <= -r)
    sel12 <- sort(union(sel1, sel2))
    
    ## post hoc bounds in selections
    n1 <- length(sel1)
    FP1 <- maxFP(pval[sel1], thr = thr)
    TP1 <- n1 - FP1
    FDP1 <- round(FP1/max(n1, 1), 2)
    
    n2 <- length(sel2)
    FP2 <- maxFP(pval[sel2], thr = thr)
    TP2 <- n2 - FP2
    FDP2 <- round(FP2/max(n2, 1), 2)
    
    n12 <- length(sel12)
    FP12 <- maxFP(pval[sel12], thr = thr)
    TP12 <- n12 - FP12
    FDP12 <- round(FP12/max(n12, 1), 2)
    
    ## graphical parameters
    cols <- rep(col[1], m)
    cols[c(sel1, sel2)] <- col[2]
    
    cexs <- rep(cex[1], m)
    cexs[sel12] <- cex[2]

    xlab <- "Fold change (log scale)"
    ylab <- bquote("p-value (-" ~ log[10] ~ "scale)")
    infty <- 100
    if (is.null(ylim)) {
        ylim <- c(0, max(logp))
    }
    plot(fc, logp, pch = pch, cex = cexs, col = cols, 
         xlab = xlab, ylab = ylab, ylim = ylim)
    # axis4 <- thrYaxis(thr, max(ylim))
    # axis(side = 4, at = axis4$pvalue, labels = axis4$num, las = 1)
    # abline(h = -log10(thr[1:100]), col = "lightgray")
    rect(xleft = -infty, ybottom = y_thr, xright = -r, ytop = infty, 
         col = col[3], border = NA, lwd = 2)
    rect(xleft = r, ybottom = y_thr, xright = infty, ytop = infty, 
         col = col[3], border = NA, lwd = 2)
    abline(h = y_thr, col = "gray")
    abline(v = c(-1, 1)*r, col = "gray")

    if (bounds) {
        bq <- bquote(atop(.(n1) ~ "genes", 
                          "TP"  >= .(TP1) ~ ";  FDP" <= .(FDP1)))
        legend("topright", legend = bq, border = "white", bty = "n", text.col = 1)
        bq <- bquote(atop(.(n2) ~ "genes",
                          "TP"  >= .(TP2) ~ "; FDP" <= .(FDP2)))
        legend("topleft", legend = bq, border = "white", bty = "n", text.col = 1)
        
        bq <- bquote(atop(.(n12) ~ "genes selected",
                          "At least" ~ .(TP12) ~ "true positives (FDP" <= .(FDP12) ~")"))
        title(bq)
    }
    invisible(sel12)
}
