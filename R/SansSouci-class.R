#' Create an object of class 'SansSouci'
#'
#' @param Y A matrix of \eqn{m} variables (hypotheses) by \eqn{n} observations
#' @param groups A numeric vector of \eqn{n} values in \eqn{0, 1}, the groups of
#'   observations on which to perform two-sample tests
#' @param truth An optional numeric vector of $m$ values in ${0,1}$, the status
#'   of each null hypothesis (0 means H0 is true, 1 means H1 is true). Typically
#'   used in simulations.
#' @return An object of class "SansSouci"
#' @export
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
#' res <- fit(a, B = 100, alpha = 0.1, family="Beta", K=10)
#' label(res)
#' 
SansSouci <- function(Y, groups, truth = NULL) {
    ugroups <- unique(groups)
    n_groups <- length(ugroups)
    
    if (n_groups > 1) {
        categCheck(groups, ncol(Y))
    }
    if (!is.null(truth)) {
        categCheck(truth, nrow(Y))
    }
    input = list(Y = Y, 
                 groups = groups, 
                 n_groups = n_groups,
                 m = nrow(Y))
    input$truth <- truth
    obj <- structure(list(input = input,
                          parameters = NULL,
                          output = NULL), 
                     class = "SansSouci")
    obj
}


#' Create an object of class 'SansSouci' from simulations
#' 
#' Create a 'SansSouci' object from simulation in the Gaussian 
#' equi-correlated model
#' 
#' @param ... Parameters to be passed to [gaussianSamples]
#' @seealso [gaussianSamples]
#' @export
#' @examples
#' obj <- SansSouciSim(m = 543, rho = 0.4, n = 210,
#'                         pi0 = 0.8, SNR = 3, prob = 0.5)
#' alpha <- 0.1
#' 
#' # Adaptive Simes (lambda-calibration)
#' res <- fit(obj, B = 100, alpha = alpha, family = "Simes")
#' res
#' # upper bound on number of signals if the entire data set
#' # (and corresponding lower bound on FDP)
#' predict(res)
#' 
#' # confidence curve 
#' plot(res)
#' 
#' # comparison to other confidence curves
#' # Parametric Simes (no calibration -- assume positive dependence (PRDS))
#' res0 <- fit(obj, B = 0, alpha = alpha, family = "Simes")
#' res0
#' 
#' # Oracle
#' oracle <- fit(obj, alpha = alpha, family = "Oracle")
#' oracle
#' 
#' confs <- list(Simes = predict(res0, all = TRUE),
#'               "Simes+calibration" = predict(res, all = TRUE),
#'               "Oracle" = predict(oracle, all = TRUE))
#' plotConfCurve(confs)
SansSouciSim <- function(...) {
    sim <- gaussianSamples(...)
    SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H)
}

#' Basic methods for class `SansSouci`
#' 
#' @param object An object of class \code{SansSouci}
#' @name SansSouci-methods
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
#' label(res)
#' print(res)
#' str(pValues(res))
#' str(foldChanges(res)) 
#' str(thresholds(res))
NULL
#> NULL

#' Generic functions for S3 class SansSouci
#' 
#' @name all-generics
#' @param object An object. See individual methods for specifics
NULL
#> NULL

#' @describeIn all-generics ' Get the number of hypotheses
#' @param object An object. See individual methods for specifics
#' @export
nHyp <- function(object) UseMethod("nHyp")

#' `nHyp`: get the number of hypotheses
#' 
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
nHyp.SansSouci <- function(object) {
    object$input$m
}

#' @describeIn all-generics Get the number of observations
#' @param object An object. See individual methods for specifics
#' @export
nObs <- function(object) UseMethod("nObs")

#' `nObs` Get the number of observations
#' 
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
nObs.SansSouci <- function(object) {
    ncol(object$input$Y)
}

#' @describeIn all-generics Get the label of hypotheses tested
#' @param object An object. See individual methods for specifics
#' @export
label <- function(object) UseMethod("label")

#' `label` Get the label of a post hoc method
#' 
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
label.SansSouci <- function(object) {
    param <- object$parameters
    if (is.null(param)) {
        return(NULL)
    }
    lab <- param$family
    if (!(lab %in% c("Simes", "Oracle")) & (!is.null(param$K))) {
        lab <- sprintf("%s(K=%d)", lab, param$K)
    } 
    return(lab)
}

#' Print 'SansSouci' objects
#' 
#' @rdname SansSouci-methods
#' @param x An object of class `SansSouci`
#' @param ... Not used
#' @param verbose Should detailed output be printed? Defaults to FALSE
#' @importFrom utils str
#' @export
print.SansSouci <- function(x, ..., verbose = FALSE) {
    object <- x; rm(x)
    cat("'SansSouci' object:\n")
    input <- object$input
    if (!is.null(input)) {
        cat("\tNumber of hypotheses: ", nHyp(object), "\n")
        cat("\t", input$n_group, "-sample data", "\n", sep="")
        cat("\n")
        if (verbose) {
            cat("Data:")
            cat("\n")
            str(input$Y)
            cat("\n")
        }
        truth <- input$truth
        if (!is.null(truth)) {
            cat("Truth: ")
            cat("\n")
            cat("\t", sum(truth), "false null hypotheses (signals)")
            pi0 <- 1 - mean(truth)
            cat(" out of ", nHyp(object), " (pi0=", round(pi0, 3), ")", sep = "")
            cat("\n")
        }
    }
    params <- object$parameters
    if (!is.null(params)) {
        cat("Parameters:", "\n")
        cat("\tTest function:", params$funName, "\n")
        cat("\tNumber of permutations: B=", params$B, "\n", sep = "")
        cat("\tSignificance level: alpha=", params$alpha, "\n", sep = "")
        cat("\tReference family:", params$fam, "\n")
        cat("\t(of size: K=", params$K, ")", "\n", sep = "")
        cat("\n")
    }
    output <- object$output
    if (!is.null(output)) {
        cat("Output:\n")
        if (verbose) {
            cat("\tp-values:\n")
            print(summary(output$p.values))
            cat("\tCalibrated thresholds:")
            if (all(output$thr %in% c(0,1))) {
                print(table(output$thr))
                cat("\n")
            } else {
                print(summary(output$thr))
                cat("\tPivotal statistic:\n")
                print(summary(output$piv_stat))
            }
            cat("\n")
        }
        cat("\tCalibration parameter: lambda=", output$lambda, "\n", sep = "")
    }
    invisible(object)
}

#' @importFrom generics fit
#' @export 
generics::fit

#' @describeIn calibrateJER Fit SansSouci object
#' @param object An object of class `SansSouci`
#' @param family A character value which can be \describe{
#'   \item{Simes}{The classical family of thresholds introduced by Simes (1986):
#'   \eqn{\alpha*k/m}. This family yields JER control if the test
#'   statistics are positively dependent (PRDS) under H0.}
#'
#'    \item{Beta}{A family of thresholds that achieves "balanced" JER control under independence}
#'    
#'   \item{Oracle}{A family such that the associated bounds correspond to the true numbers/proportions of true/false positives. "truth" must be available in object$input$truth.}
#'   }
#' @param ... Not used
#' @return A 'fitted' object of class 'SansSouci'. It is a list of three elements
#'  - input: see [SansSouci]
#'  - param: the input parameters, given as a list
#'  - output: the outputs of the calibration, see [calibrateJER]
#' @export
fit.SansSouci <- function(object, alpha, B = ceiling(10/alpha),
                          alternative = c("two.sided", "less", "greater"),
                          rowTestFUN = NULL, 
                          family = c("Simes", "Beta", "Oracle"), 
                          maxStepsDown = 10L, K = nHyp(object), 
                          flavor = c("v0", "v1"),
                          verbose = TRUE, ...) {
    alternative <- match.arg(alternative)
    family <- match.arg(family)
    flavor <- match.arg(flavor)
    if (family == "Oracle") {
        truth <- object$input$truth
        if (is.null(truth)) {
            stop("'truth' should be available for 'Oracle'. See ?SansSouci")
        }
    }
    Y <- object$input$Y
    groups <- object$input$groups
    n_groups <- object$input$n_groups
    m <- nHyp(object)
    funName <- NA_character_
    if (is.null(rowTestFUN)) {
        if (n_groups == 1) {
            rowTestFUN <- function(mat, categ, alternative) {
                T <- rowSums(mat)/sqrt(length(categ))
                p <- switch(alternative, 
                            "two.sided" = 2*(pnorm(abs(T), lower.tail = FALSE)),
                            "greater" = pnorm(T, lower.tail = FALSE),
                            "less" = pnorm(T, lower.tail = TRUE))
                data.frame(statistic = T, parameter = NA, p.value = p)
            }
            funName <- "testBySignFlipping"
        } else if (n_groups == 2) {
            rowTestFUN <- rowWelchTests
            funName <- "rowWelchTests"
        }
    }  else {
        funName <- as.character(substitute(rowTestFUN))
    }
    object$parameters <- list(
        alpha = alpha,
        B = B,
        alternative = alternative, 
        rowTestFUN = rowTestFUN, 
        funName = funName,
        family = family, 
        maxStepsDown = maxStepsDown, 
        K = K)
    
    if (flavor == "v0") {
        if (B > 0 && family != "Oracle") {
            cal <- calibrateJER(X = Y, categ = groups, B = B, alpha = alpha, 
                                alternative = alternative, 
                                rowTestFUN = rowTestFUN, 
                                refFamily = family, 
                                maxStepsDown = maxStepsDown, 
                                K = K, verbose = verbose)
        } 
        else {
            # no calibration!
            rwt <- rowTestFUN(mat = Y, categ = groups, alternative = alternative)
            p.values <- rwt$p.value
            fc <- rwt$meanDiff
            lambda <- alpha
            if (family %in% c("Simes", "Linear")) {
                thr <- SimesThresholdFamily(m, kMax = K)(alpha)
            } else if (family == "Beta") {
                thr <- BetaThresholdFamily(m, kMax = K)(alpha)
            } else if (family == "Oracle") {
                lambda <- 0  # 100% confidence -- should it be set to alpha for consistency?
                thr <- object$input$truth
            }
            cal <- list(p.values = p.values,
                        fold_changes = fc,
                        piv_stat = NULL,
                        thr = thr,
                        lambda = alpha,
                        conf_bound = NULL)
        }
    } else if (flavor == "v1") {
        rwt <- rowWelchTests(X, groups)
        p.values <- rwt$p.value
        cal <- list(p.values = p.values,
                    fold_changes = rwt$meanDiff)
        if (B > 0 && family != "Oracle") {
            p0 <- sansSouci:::get_perm_p2(X = X, categ = groups, 
                                          B = B, rowTestFUN = rowTestFUN)
            
            calib <- get_calibrated_thresholds_sd(p0 = p0, m = m, alpha = alpha, 
                                                  family = family, K = K, 
                                                  p = p.values, maxStepsDown = maxStepsDown)
            cal$piv_stat <- calib$pivStat
            cal$thr <- calib$thr
            cal$lambda <- calib$lambda
        } else {
            lambda <- alpha
            if (family == "Linear") {
                thr <- t_linear(alpha, seq_len(m), m)
            } else if (family == "Beta") {
                thr <- t_beta(alpha, seq_len(m), m)
            } else if (family == "Oracle") {
                lambda <- 0  # 100% confidence -- should it be set to alpha for consistency?
                thr <- object$input$truth
            }
            cal$thr <- cal
        }
    }
    object$output <- cal
    object
}

#' @describeIn all-generics Get p-values
#' @param object An object. See individual methods for specifics
#' @export
pValues <- function(object) UseMethod("pValues")

#' `pValues`: get p-values
#' 
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
pValues.SansSouci <- function(object) {
    object$output$p.values
}

#' @describeIn all-generics Get fold changes
#' @param object An object. See individual methods for specifics
#' @export
foldChanges <- function(object) UseMethod("foldChanges")

#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
foldChanges.SansSouci <- function(object) {
    object$output$fold_changes
}

#' @describeIn all-generics Get thresholds
#' @param object An object. See individual methods for specifics
#' @export
thresholds <- function(object) UseMethod("thresholds")

#' `thresholds`: get thresholds
#' 
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
thresholds.SansSouci <- function(object) object$output$thr

#' Plot confidence bound on the true/false positives among most significant items
#' 
#' @param x An object of class 'SansSouci'
#' @param y Not used
#' @param xmax Right limit of the plot
#' @param ... Further arguments to be passed to \code{bound}
#' @export
plot.SansSouci <- function(x, y, xmax = nHyp(x), ...) {
    cb <- predict(x, all = TRUE, ...)
    plotConfCurve(cb, xmax = xmax)
}

