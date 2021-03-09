#' Create an object of class "SansSouci"
#' 
#' @param Y A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param groups An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#' @return An object of class "SansSouci"
#' @export
#' 
#' @examples
#' data(expr_ALL, package = "sansSouci.data")
#' groups <- ifelse(colnames(expr_ALL)=="NEG", 0, 1)
#' table(groups)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#' print(a)
#' nHyp(a)
#' nObs(a)
#' 
#' res <- fit(a, B = 100, alpha = 0.1)
#' print(res)
#' label(res)
#' 
#' res <- fit(a, B = 100, alpha = 0.1, refFamily="Beta", K=10)
#' label(res)
SansSouci <- function(Y, groups) {
    obj <- structure(list(input = list(Y = Y, 
                                       groups = groups),
                          parameters = NULL,
                          alpha = NULL,
                          output = NULL), 
                     class = "SansSouci")
    return(obj)
}

#' SansSouci class
#'
#' @param object An object of class \code{SansSouci}
#' @name SansSouci-class
#' @examples
#' data(expr_ALL, package = "sansSouci.data")
#' groups <- ifelse(colnames(expr_ALL)=="NEG", 0, 1)
#' table(groups)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#' print(a)
#' nHyp(a)
#' nObs(a)
#' label(a)
#' 
#' res <- fit(a, B = 100, alpha = 0.1)
#' print(res)
#' label(res)
NULL
#> NULL


#' Get the number of hypotheses
#' 
#' @param object An object. See individual methods for specifics
#' @param ... Other arguments passed to methods
#' @export
nHyp <- function(object) UseMethod("nHyp")

#' @rdname SansSouci-class
#' @export
nHyp.SansSouci <- function(object) {
    nrow(object$input$Y)
}

#' Get the number of observations
#' 
#' @inheritParams nHyp
#' @export
nObs <- function(object) UseMethod("nObs")

#' @rdname SansSouci-class
#' @export
nObs.SansSouci <- function(object) {
    ncol(object$input$Y)
}

#' Get the label of hypotheses tested
#' 
#' @inheritParams nHyp
#' @export
label <- function(object) UseMethod("label")

#' @rdname SansSouci-class
#' @export
label.SansSouci <- function(object) {
    param <- object$parameters
    if (is.null(param)) {
        return(NULL)
    }
    lab <- param$refFamily
    if (param$K < nHyp(object)) {
        lab <- sprintf("%s(K=%d)", lab, param$K)
    } 
    return(lab)
}

#' @rdname SansSouci-class
#' @export
print.SansSouci <- function(x, ...) {
    object <- x; rm(x)
    cat("'SansSouci' object")
    input <- object$input
    if (!is.null(input)) {
        msg <- sprintf(" with %d hypotheses.", 
                       nHyp(object))
        cat(msg)
        cat("\n")
        cat("Data:")
        cat("\n")
        str(input$Y)
        cat("\n")
    }
    params <- object$parameters
    if (!is.null(params)) {
        cat("Parameters:", "\n")
        cat("\tTest function:", params$funName, "\n")
        cat("\tNumber of permutations: B=", params$B, "\n", sep = "")
        cat("\tSignificance level: alpha=", params$alpha, "\n", sep = "")
        cat("\tReference family:", params$refFam, "\n")
        cat("\t(of size: K=", params$K, ")", "\n", sep = "")
        cat("\n")
    }
    output <- object$output
    if (!is.null(output)) {
        cat("Output:\n")
        cat("\tp-values:\n")
        print(summary(output$p.values))
        cat("\tPivotal statistic:\n")
        print(summary(output$pivStat))
        cat("\tCalibration parameter: lambda=", output$lambda, "\n", sep = "")
        cat("\tCalibrated thresholds:\n")
        print(summary(output$thr))
        cat("\n")
    }
    invisible(object)
}

#' @describeIn calibrateJER Fit SansSouci object
#' @importFrom generics fit
#' @export fit
#' @export
fit.SansSouci <- function(object, B, alpha, 
                           alternative = c("two.sided", "less", "greater"),
                           rowTestFUN = rowWelchTests, 
                           refFamily = c("Simes", "Beta"), 
                           maxStepsDown = 10L, K = nHyp(object), 
                           verbose = TRUE, ...) {
    alternative <- match.arg(alternative)
    refFamily <- match.arg(refFamily)
    object$parameters <- list(B = B,
                     alpha = alpha,
                     alternative = alternative, 
                     rowTestFUN = rowTestFUN, 
                     funName = as.character(substitute(rowTestFUN)),
                     refFamily = refFamily, 
                     maxStepsDown = maxStepsDown, 
                     K = nHyp(object), verbose = verbose)
    Y <- object$input$Y
    groups <- object$input$groups
    
    if (B>0) {
        cal <- calibrateJER(X = Y, categ = groups, B=B, alpha=alpha, 
                        alternative = alternative, 
                        rowTestFUN = rowTestFUN, 
                        refFamily = refFamily, 
                        maxStepsDown = maxStepsDown, 
                        K = K, verbose = verbose)
    } else {
        ## B=0: no calibration!
        rwt <- rowTestFUN(mat = Y, categ = groups, alternative = alternative)
        p.values <- rwt$p.value
        fc <- rwt$meanDiff
        lambda <- alpha
        if (refFamily %in% c("Simes", "Linear")) {
            thr <- SimesThresholdFamily(m, kMax = K)(alpha)
        } else if (refFamily == "Beta") {
            thr <- BetaThresholdFamily(m, kMax = K)(alpha)
        }
        cal <- list(p.values = p.values,
                    fold_changes = fc,
                    piv_stat = NULL,
                    thr = thr,
                    lambda = alpha,
                    conf_env = NULL)
    }
    object$output <- cal
    object
}


#' Get p-values
#' 
#' @inheritParams nHyp
#' @export
pValues <- function(object, ...) UseMethod("pValues")

#' @rdname SansSouci-class
#' @export
pValues.SansSouci <- function(object, ...) {
    object$output$p.values
}

#' Get fold changes
#' 
#' @inheritParams nHyp
#' @export
fold_changes <- function(object, ...) UseMethod("fold_changes")

#' @rdname SansSouci-class
#' @export
fold_changes.SansSouci <- function(object, ...) object$output$fold_changes

#' Get thresholds
#' @export
#' @inheritParams nHyp
thresholds <- function(object, ...) UseMethod("thresholds")

#' @rdname SansSouci-class
#' @export
thresholds.SansSouci <- function(object, ...) object$output$thr

#' Plot confidence envelope on the true/false positives among most significant items
#' 
#' @param x An object of class 'SansSouci'
#' @param y Not used
#' @param ... Further arguments to be passed to \code{bound}
#' @export
plot.SansSouci <- function(x, y, xmax = nHyp(x), ...) {
    ce <- bound(x, envelope = TRUE, ...)
    plotConfEnvelope(ce, xmax = xmax)
}

