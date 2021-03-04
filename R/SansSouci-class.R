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
SansSouci <- function(Y, groups) {
    obj <- structure(list(input = list(Y = Y, 
                                       groups = groups),
                          parameters = NULL,
                          alpha = NULL,
                          output = NULL), 
                     class = "SansSouci")
    obj
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
    }
    params <- object$parameters
    if (!is.null(params)) {
        cat("Parameters:")
        cat("\n")
        str(params)
    }
    output <- object$output
    if (!is.null(output)) {
        cat("Output:")
        cat("\n")
        str(output)
    }
    invisible(object)
}

#' @describeIn calibrateJER Fit SansSouci object
#' @importFrom generics fit
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
                     refFamily = refFamily, 
                     maxStepsDown = maxStepsDown, 
                     K = K, verbose = verbose)
    Y <- object$input$Y
    groups <- object$input$groups
    
    cal <- calibrateJER(X = Y, categ = groups, B=B, alpha=alpha, 
                        alternative = alternative, 
                        rowTestFUN = rowTestFUN, 
                        refFamily = refFamily, 
                        maxStepsDown = maxStepsDown, 
                        K = K, verbose = verbose)
    object$output <- cal
    object
}


#' @export
p_values <- function(object) UseMethod("p_values")
#' @export
p_values.SansSouci <- function(object) object$output$p.values

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
