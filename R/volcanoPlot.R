#' Volcano plot
#' 
#' @param x An object. See individual methods for specifics
#' @param ... Other arguments passed to methods
#' @export
volcanoPlot <- function(x, ...) UseMethod("volcanoPlot")

#' @inheritParams volcanoPlot
#' @export
volcanoPlot.SansSouci <- function(x, ...) {
    object <- x; rm(x);
    pval <- p_values(object)
    fc <- fold_changes(object)
    thr <- thresholds(object)
    volcanoPlot(pval = pval, fc = fc, thr = thr, ...)
}

#' Volcano plot
#' 
#' Volcano plot for differential expression studies
#' 
#' @param X A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations

#' @param categ A numeric vector of \code{n} categories in \eqn{0, 1} for the
#'   observations

#' @param thr A numeric vector of length K, a JER controlling family
#' 
#' @param p A numeric value, the p-value threshold under which genes are selected
#' @param q A numeric value, the q-value (or FDR-adjusted p-value) threshold under which genes are selected
#' @param r A numeric value, the absolute fold change above which genes are selected
#' @param cex A numeric vector of length 2, the relative magnification factor for unselected (\code{cex[1]}) and unselected (\code{cex[2]}) genes. 
#' 
#' @param col A vector of length 3
#' @param pch An integer or single character string specifying the plotting character, see \code{\link{par}}
#' @param ylim A numeric vector of length 2, the \eqn{y} limits of the plot
#'
#' @details A Welch T-test of differential expression between the two categories
#'   defined by \code{categ} are applied for each gene using the
#'   \code{\link{rowWelchTests}} function, which also outputs the "fold change"
#'   (mean difference in log scale) between the two categories.
#' @return The indices of selected genes (returned invisibly)
#' 
#' @export
#' @importFrom graphics abline legend rect title
#' @importFrom stats p.adjust
#' @seealso Volcano plot shiny app at \url{https://pneuvial.shinyapps.io/volcano-plot/}
#' 
#' @examples
#' m <- 500
#' pi0 <- 0.5
#' m1 <- m-m*pi0
#' SNR <- 5*(runif(m1)-0.5)
#' sim <- gaussianSamples(m = m, rho = 0.2, n = 100,
#'                        pi0 = pi0, SNR = SNR, prob = 0.5)
#' X <- sim$X
#' categ <- sim$categ
#' alpha <- 0.2
#' cal <- calibrateJER(X = X, categ = categ, B = 1e2, alpha = alpha, refFamily="Simes")
#' tests <- rowWelchTests(X, categ)
#' sel <- volcanoPlot(tests$p.value, tests$meanDiff, thr = cal$thr, q = 0.2, r = 0.2, ylim = c(0, 6))
#' 
#' # Compare bound to reality
#' TP <- sum(sim$H[sel])
#' FP <- sum(1-sim$H[sel])
#' FDP <- FP/length(sel)

volcanoPlot.numeric <- function(pval, fc, thr,
                        p = 1, q = 1, r = 0,
                        cex = c(0.2, 0.6), 
                        col = c("#33333333", "#FF0000", "#FF666633"),
                        pch = 19, ylim = NULL, ...) {
    if (p <1 && q < 1) {
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
        y_thr <- min(logp[y_sel])              ## threshold on the log(p-value) scale
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
    bq <- bquote(atop(.(n1) ~ "genes", 
                      "TP"  >= .(TP1) ~ ";  FDP" <= .(FDP1)))
    legend("topright", legend = bq, border = "white", bty = "n", text.col = 1)
    bq <- bquote(atop(.(n2) ~ "genes", 
                       "TP"  >= .(TP2) ~ "; FDP" <= .(FDP2)))
    legend("topleft", legend = bq, border = "white", bty = "n", text.col = 1)
    
    bq <- bquote(atop(.(n12) ~ "genes selected", 
                      "At least" ~ .(TP12) ~ "true positives (FDP" <= .(FDP12) ~")"))
    title(bq)
    invisible(sel12)
}

thrYaxis <- function(thr, maxlogp){
    df1 <- data.frame(num = 1:length(thr)-1, pvalue=-log10(thr))
    df2 <- data.frame(df1[c(1),])
    valeurTest <- df2[c(dim(df2)[1]),"pvalue"]
    for (i in 1:dim(df1)[1]){
        mod <- if(df1[i,"num"] < 100){ 1} else if(df1[i,"num"] < 500){ 10} else if(df1[i,"num"] < 1000){50}else{100}
        if (valeurTest - df1[i,"pvalue"] > 0.3*maxlogp/12.5 & df1[i,"num"]%%mod == 0){
            df2 <- rbind(df2, (df1[i,]))
            valeurTest <- df1[i,"pvalue"]
        }
    }
    return(df2)
}
