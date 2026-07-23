#' Volcano plot
#'
#' Volcano plot for differential expression studies
#'
#' @rdname volcanoPlot
#' @param x An object of class `SansSouci`, or alternatively a numeric vector
#'   of fold changes when calling `volcanoPlot.numeric` directly (x axis of
#'   the volcano plot). See [volcanoPlot.numeric()] or [volcanoPlot.SansSouci()]
#'   for specifics
#' @param ... Other arguments passed to methods
#' @export
volcanoPlot <- function(x, ...) UseMethod("volcanoPlot")


#' Volcano plot for a numeric vector of fold changes and p-values
#' @param x A numeric vector of fold changes (x axis of the volcano plot)
#' @param p_value A numeric vector of p-values, of the same length as `x`
#'   (y axis of the volcano plot)
#' @param thr A numeric vector of length K, a JER controlling family, used to estimate post-hoc bounds
#' @param p_value_bound A numeric vector of p-values, of the same length as `x`, used to estimate post-hoc bounds. Defaults to 'p_value'
#' @param p A numeric value, the p-value threshold under which features are selected
#' @param q A numeric value, the q-value (or FDR-adjusted p-value) threshold under which features are selected
#' @param r A numeric value, the absolute fold change above which features are selected
#' @param cex A numeric vector of length 2, the relative magnification factor for unselected (\code{cex[1]}) and unselected (\code{cex[2]}) features.
#'
#' @param col A vector of length 3
#' @param pch An integer or single character string specifying the plotting character, see \code{\link{par}}
#' @param feature_label A character, the label to be used to designate individual hypotheses in the plot title
#' (e.g. "gene", "proteins", ...)
#' @param ylim A numeric vector of length 2, the \eqn{y} limits of the plot
#' @param add_signed_selections A boolean value: should the post hoc bounds for the subselections corresponding to positive and negative fold change be displayed? Defaults to TRUE
#' @param ... Not used
#' @details The p-values play two distinct roles here: they are displayed as the
#'   y axis of the volcano plot, and they are also used to compute post hoc
#'   bounds. In general the same p-values are used. See the vignette
#'   <https://sanssouci-org.github.io/sanssouci/articles/post-hoc_differential-expression_RNAseq.html#custom-statistics-example-using-limma-voom>
#'   for an example where two different sets of p-values are used.
#' @return A 'ggplot' object containing the volcano plot, with the indices of selected features returned as an attribute named 'selection'
#'
#' @exportS3Method
#' @importFrom graphics abline legend rect title
#' @importFrom stats p.adjust
#' @export
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' fold_changes <- foldChanges(res)
#' p_values <- pValues(res)
#' thr <- thresholds(res)
#' volcanoPlot(x = fold_changes[1, ], p_value = p_values[1, ], thr = thr, 
#'   q = 0.2, r = 0.2, ylim = c(0, 4), feature_label = "gene")
#'
#' @seealso Volcano plot shiny app at \url{ https://shiny-iidea-sanssouci.apps.math.cnrs.fr/}
#' @seealso Volcano plot for objects of class 'SansSouci': [volcanoPlot.SansSouci()]
volcanoPlot.numeric <- function(x, p_value, thr, p_value_bound = p_value, 
                                p = 1, q = 1, r = 0,
                                cex = c(0.4, 1.5),
                                col = c("#33333333", "#FF0000", "#FF666633"),
                                pch = 19, feature_label = "feature",
                                ylim = NULL, add_signed_selections = TRUE,
                                ...) {
  fold_change <- x
  if (p < 1 && q < 1) {
    warning("Filtering both on p-values and BH-adjusted p-values")
  }
  m <- length(p_value)

  ## sanity checks
  stopifnot(length(fold_change) == m)
  stopifnot(length(p_value) == m)
  stopifnot(length(p_value_bound) == m)
  stopifnot(length(thr) <= m)
  
  logp <- -log10(p_value)
  adjp <- p.adjust(p_value, method = "BH") ## adjusted p-values
  y_sel <- which((adjp <= q) &       ## selected by q-value
                   (p_value <= p))         ## and/or p-value
  y_thr <- Inf
  if (length(y_sel) > 0) {
    y_thr <- min(logp[y_sel]) ## threshold on the log(p-value) scale
  }

  ## feature selections
  sel1 <- which(logp >= y_thr & fold_change >= r)
  sel2 <- which(logp >= y_thr & fold_change <= -r)
  sel12 <- sort(union(sel1, sel2))

  ## post hoc bounds in selections
  n1 <- length(sel1)
  FP1 <- maxFP(p_value_bound[sel1], thr = thr)
  TP1 <- n1 - FP1
  FDP1 <- round(FP1 / max(n1, 1), 2)

  n2 <- length(sel2)
  FP2 <- maxFP(p_value_bound[sel2], thr = thr)
  TP2 <- n2 - FP2
  FDP2 <- round(FP2 / max(n2, 1), 2)

  n12 <- length(sel12)
  FP12 <- maxFP(p_value_bound[sel12], thr = thr)
  TP12 <- n12 - FP12
  FDP12 <- round(FP12 / max(n12, 1), 2)

  ## graphical parameters
  cols <- rep(col[1], m)
  cols[c(sel1, sel2)] <- col[2]

  cexs <- rep(cex[1], m)
  cexs[sel12] <- cex[2]

  xlab <- "Fold change (log scale)"
  ylab <- bquote("p-value (-" ~ log[10] ~ "scale)")


  df <- data.frame(
    log_pval = logp, logfc = fold_change,
    selected = factor(ifelse(seq(m) %in% sel12, "In", "Out"),
                      levels = c("Out", "In")
    )
  )
  title <- sprintf(
    "%d %s selected\nAt least %d true positives (FDP \u2264 %.2f)",
    n12, pluralize(word = feature_label, n = n12), TP12, FDP12
  )

  vp <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$logfc, y = .data$log_pval, color = .data$selected)
  ) +
    ggplot2::geom_hline(
      yintercept = y_thr,
      linetype = "dashed",
      color = "darkgrey"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1) * r,
      linetype = "dashed",
      color = "darkgrey"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      alpha = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "solid",
      color = "black",
      alpha = 0.4
    ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$selected)) +
    ggplot2::scale_color_manual(values = col[1:2]) +
    ggplot2::scale_size_manual(values = cex) + # cex = c(0.2, 0.6)
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    # color selected zone
    ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = -r,
      ymin = y_thr, ymax = Inf,
      fill = col[3], alpha = 0.3
    ) +
    ggplot2::annotate(
      "rect",
      xmin = r, xmax = Inf,
      ymin = y_thr, ymax = Inf,
      fill = col[3], alpha = 0.3
    )

  if (add_signed_selections) {
    txt_right <- sprintf(
      "%d %s\nTP \u2265 %d ; FDP \u2264 %.2f",
      n1, pluralize(word = feature_label, n = n1), TP1, FDP1
    )

    txt_left <- sprintf(
      "%d %s\nTP \u2265 %d ; FDP \u2264 %.2f",
      n2, pluralize(word = feature_label, n = n2), TP2, FDP2
    )
    # bounds
    vp <- vp +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        label = txt_right,
        hjust = 1.05,
        vjust = 1.2
      ) +
      ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        label = txt_left,
        hjust = -0.05,
        vjust = 1.2
      )
  }
  attr(vp, "selection") <- sel12
  return(vp)
}


#' Volcano plot for a `SansSouci` object
#' 
#' @param fold_change An optional vector of fold changes, of the same length as `nHyp(x)`, used for volcanoPlot x-axis. If not specified, `foldChanges(x)` is used.
#' @param p_value A vector of p-values, of the same length as `nHyp(x)`, used for volcanoPlot y-axis. If not specified, `pValues(x)` is used
#' @param contrast_name A character value, the selected contrast. Should be chosen in `x$input$contrast_name`.
#' @inheritParams volcanoPlot.numeric
#' @details The default is to use the fold changes and p-values from the input SansSouci object. See the vignette
#'   <https://sanssouci-org.github.io/sanssouci/articles/post-hoc_differential-expression_RNAseq.html#custom-statistics-example-using-limma-voom>
#'   for an example where custom fold changes and p-values are used.
#' 
#' @export
#'
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' volcanoPlot(res, q = 0.2, r = 0.2, ylim = c(0, 4), feature_label = "gene")
volcanoPlot.SansSouci <- function(x,
                                  fold_change = foldChanges(x)[contrast_name, ],
                                  p_value = pValues(x)[contrast_name, ],
                                  p = 1, q = 1, r = 0,
                                  contrast_name = x$input$contrast_name[1],
                                  cex = c(0.4, 1.5),
                                  col = c("#33333333", "#FF0000", "#FF666633"),
                                  pch = 19, feature_label = "feature",
                                  ylim = NULL, ...) {
  object <- x
  if (!(contrast_name %in% rownames(pValues(object)))) {
    stop(paste(
      "Choose a contrast in",
      paste(rownames(pValues(object)),
            collapse = ", "
      )
    ))
  }
  fold_change <- force(fold_change)
  p_value <- force(p_value)
  
  if (object$input$type == "1 sample") {
    stop("Can't do a volcano plot for one-sample tests!")
  }
  m <- object$input$n_dimensions
  stopifnot(m == length(fold_change))
  stopifnot(m == length(p_value))
  p_value_bound <- pValues(object)[contrast_name, ]
  thr <- thresholds(object)[1:m] # we select at most m hypotheses here
  
  volcanoPlot(
    x = fold_change, p_value = p_value, 
    thr = thr, p_value_bound = p_value_bound,
    p = p, q = q, r = r,
    cex = cex,
    col = col,
    pch = pch,
    feature_label = feature_label,
    ylim = ylim, ...
  )
}
