#' @export
#' @importFrom graphics abline legend rect title
#' @importFrom stats p.adjust
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