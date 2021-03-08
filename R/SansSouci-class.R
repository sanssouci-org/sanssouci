#' Create object of class "SansSouci"
#' 
#' @param Y A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param groups An optional numeric vector of \eqn{n} values in \eqn{0, 1}
#' @return An object of class "SansSouci"
#' @export
#' @examples
#' data(expr_ALL, package = "sansSouci.data")
#' groups <- ifelse(colnames(expr_ALL)=="NEG", 0, 1)
#' table(groups)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#' print(a)
#' n_hyp(a)
#' n_obs(a)
#' res <- fit(a, B = 100, alpha = 0.1)
#' print(res)
#' label(res)
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

#' @export
n_hyp <- function(object) UseMethod("n_hyp")
#' @export
n_hyp.SansSouci <- function(object) nrow(object$input$Y)

#' @export
n_obs <- function(object) UseMethod("n_obs")
#' @export
n_obs.SansSouci <- function(object) ncol(object$input$Y)

#' @export
label <- function(object) UseMethod("label")
#' @export
label.SansSouci <- function(object) {
    param <- object$parameters
    lab <- param$refFamily
    if (param$K < n_hyp(object)) {
        lab <- sprintf("%s(K=%d)", lab, param$K)
    } 
    return(lab)
}

#' @export
print.SansSouci <- function(object) {
    cat("'SansSouci' object")
    input <- object$input
    if (!is.null(input)) {
        msg <- sprintf(" with %d hypotheses.", 
                       n_hyp(object))
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
                           maxStepsDown = 10L, K = n_hyp(object), 
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
                     K = n_hyp(object), verbose = verbose)
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
        p.values <- rwt$p.values
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


#' @export
p_values <- function(object) UseMethod("p_values")
#' @export
p_values.SansSouci <- function(object) object$output$p.values

#' @export
fold_changes <- function(object) UseMethod("fold_changes")
#' @export
fold_changes.SansSouci <- function(object) object$output$fold_changes

#' @export
thresholds <- function(object) UseMethod("thresholds")
#' @export
thresholds.SansSouci <- function(object) object$output$thr

#' Post hoc bounds
#' 
#' get_max_FP: upper bound on the number of false discoveries
#' 
#' @rdname posthoc-bounds
#' @param object An object of class 'SansSouci'
#' @param S A subset of indices
#' @param envelope A boolean value. Should the entire confidence envelope be returned? Defaults to FALSE. See Details.
#' @return The value of the post hoc bound on S. 
#' @details If \code{envelope} is FALSE (the default) a single value is returned, corresponding to the post hoc bound on S is returned. If \code{envelope} is TRUE, a vector of the same length as \code{S} is returned, whose i-th element is the post hoc bound on the set of i most significants items in S (in the order of \code{p_values(object)})
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020). Post hoc confidence bounds on false positives using reference families. Annals of Statistics, 48(3), 1281-1303.
#' @export
#' @examples
#'
#' m <- 543
#' pi0 <- 0.8
#' sim <- gaussianSamples(m = m, rho = 0.4, n = 100,
#'                        pi0 = pi0, SNR = 3, prob = 0.5)
#' obj <- SansSouci(Y = sim$X, groups = sim$categ)
#' obj <- fit(obj, B = 1e2, alpha = 0.1, family = "Simes")
#' 
#' S <- which(p_values(obj) < 0.02)
#' 
#' max_FP(obj, S)
#' max_FDP(obj, S)
#' min_TP(obj,S)
#' min_TDP(obj, S)
#' 
#' plot(max_FDP(obj, envelope = TRUE))
#' plot(max_FDP(obj, S, envelope = TRUE))
#' @export
max_FP <- function(object, S, envelope = FALSE, ...) UseMethod("max_FP")
#' @export
max_FP.SansSouci <- function(object, S = 1:n_hyp(object), envelope = FALSE, ...) {
    if (is.null(object$output)) {
        stop("Object not fitted yet! Please use the 'fit' function first.")
    }
    p.values <- p_values(object)
    thr <- thresholds(object)
    s <- length(S)
    sps <- sort(p.values[S])
    res <- curveMaxFP(sps, thr) ## OK to do 'thr[length(S)]' here?
    if (!envelope) {
        res <- res[s]   
    }
    res
}


#' @rdname posthoc-bounds
#' @description get_max_FDP: upper bound on the false discovery proportion
#' @export
max_FDP <- function(object, S, envelope, ...) UseMethod("max_FDP")
#' @export
max_FDP.SansSouci <- function(object, S = 1:n_hyp(object), envelope = FALSE, ...) {
    res <- max_FP(object, S = S, envelope = TRUE)/seq(along = S)
    if (!envelope) {
        res <- res[length(S)]
    }
    res <- rev(cummin(rev(res)))
    res
}


#' @rdname posthoc-bounds
#' @description get_min_TP: lower bound on the number of true discoveries
#' 
#' @export
min_TP <- function(object, S, envelope, ...) UseMethod("min_TP")
#' @export
min_TP.SansSouci <- function(object, S = 1:n_hyp(object), envelope = FALSE, ...) {
    res <- seq(along = S) - max_FP(object, S = S, envelope = TRUE)
    if (!envelope) {
        res <- res[length(S)]
    }
    res
}

#' @rdname posthoc-bounds
#' @description get_min_TDP: lower bound on the true discovery proportion
#' @export
min_TDP <- function(object, S, envelope, ...) UseMethod("min_TDP")
#' @export
min_TDP.SansSouci <- function(object, S = 1:n_hyp(object), envelope = FALSE, ...) {
    1 - max_FP(object, S = S, envelope = envelope, ...)/length(S)
}


#' Confidence envelope on the true/false positives among most significant items
#' 
#' @param object An object of class 'SansSouci'
#' @param S A subset of indices
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#' 
#' \describe{
#' \item{FP}{Upper bound on the number of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the number of true positives in the 'x' most significant items}
#' \item{FDP}{Upper bound on the proportion of false positives in the 'x' most significant items}
#' \item{TP}{Lower bound on the proportion of true positives in the 'x' most significant items}}
#' Defaults to \code{c("TP", "FDP")}.
#' 
#' @return If \code{light} is \code{FALSE} (the default), a \code{data.frame} with \eqn{m} rows and 4 columns: \describe{
#' \item{x}{Number of most significant items selected}
#' \item{label}{Label for the procedure, typically of the form '<refFamily>(<param>)'}
#' \item{bound}{Value of the post hoc bound}
#' \item{stat}{Type of post hoc bound, as specified by argument \code{bound}}. If \code{light} is \code{TRUE}, only the value of the bound is returned. 
#' }
#' @export
#' @examples
#' 
#' # Generate Gaussian data and perform multiple tests
#' sim <- gaussianSamples(m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3, prob = 0.5)
#' obj <- SansSouci(Y = sim$X, groups = sim$categ)
#' res <- fit(obj, B = 100, alpha = 0.1)
#' # calculate, print, and plot confidence envelope
#' ce <- bound(res, envelope = TRUE)
#' head(ce)
#' 
#' ce <- bound(res, S=1:10, envelope = TRUE)
#' head(ce)
#' 
#' plot(res, S=which(sim$H==1))
#' 
#' @export
bound <- function(object, S, what, envelope = FALSE, ...) UseMethod("bound")

#' @export
bound.SansSouci <- function(object, S = 1:n_hyp(object), 
                          what = c("TP", "FDP"), light = FALSE, envelope = FALSE) {
    p.values <- p_values(object)
    thr <- thresholds(object)
    lab <- label(object)
    if (max(S) > n_hyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    bounds <- conf_env(p.values = p.values[S], thr = thr, lab = lab, what = what, envelope = envelope)
    if (light) {
        bounds <- bounds[, "bound"]
    }
    return(bounds)
}


#' Plot confidence envelope on the true/false positives among most significant items
#' 
#' @param x An object of class 'SansSouci'
#' @param y Not used
#' @param ... Further arguments to be passed to \code{posthoc_bound}
#' @export
plot.SansSouci <- function(x, y, xmax = n_hyp(x), ...) {
    ce <- bound(x, envelope = TRUE, ...)
    plotConfidenceEnvelope(ce, xmax = xmax)
}

## TODO: method volcano plot !
#' @export
volcano_plot <- function(x, ...) UseMethod("volcano_plot")
#' @export
volcano_plot.SansSouci <- function(x, ...) {
    object <- x; rm(x);
    pval <- p_values(object)
    fc <- fold_changes(object)
    thr <- thresholds(object)
    volcano_plot(pval = pval, fc = fc, thr = thr, ...)
}
