#' Create an object of class 'SansSouciStruct'
#' 
#' 
#' @param struct A (dyadic) structure
#' @param leaves A list of leaves
#' @param truth An optional numeric vector of $m$ values in ${0,1}$, the status of each null hypothesis (0 means H0 is true, 1 means H1 is true). Typically used in simulations.
#' @export
#' @examples
#' s <- 100
#' q <- 7
#' m <- s*2^q
#' 
#' dd <- dyadic.from.window.size(m, s, method = 2)
#' obj1 <- SansSouciStruct(dd$C, dd$leaf_list)
#' obj <- SansSouciDyadic(m, s, flavor = "window.size", direction = "top-down")
SansSouciStruct <- function(struct, leaves, truth = NULL) {
    m <- length(unlist(leaves))
    input = list(struct = struct, 
                 leaves = leaves,
                 m = m)
    input$truth <- truth
    obj <- structure(list(input = input,
                          parameters = NULL,
                          alpha = NULL,
                          output = NULL), 
                     class = "SansSouciStruct")
    obj
}

#' @export
SansSouciDyadic <- function(m, leafSize = NULL, maxHeight = NULL,
                            flavor = c("window.size", "height", "max.height"), 
                            direction = c("bottom-up", "top-down"),
                            ...) {
    flavor <- match.arg(flavor)
    direction <- match.arg(direction)
    method <- switch(direction, 
                     "bottom-up" = 1,
                     "top-down" = 2)
    dd <- NULL
    if (flavor == "window.size") {
        dd <- dyadic.from.window.size(m, leafSize, method)
    } else if (flavor == "height") {
        dd <- dyadic.from.height(m, maxHeight, method)
    } else if (flavor == "max.height") {
        dd <- dyadic.from.max.height(m, method)
    }
    SansSouciStruct(struct = dd$C, leaves = dd$leaf_list, ...)
}


#' `nHyp`: get the number of hypotheses
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
nHyp.SansSouciStruct <- function(object) {
    object$input$m
}

#' `nLeaves`: get the number of leaves
#' 
#' @rdname all-generics
#' @param object An object. See individual methods for specifics
#' @export
nLeaves <- function(object) UseMethod("nLeaves")

#' `nLeaves`: get the number of leaves
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
nLeaves.SansSouciStruct <- function(object) {
    length(object$input$leaves)
}


#' Fit SansSouciStruct object
#' 
#' @export
#' @examples 
#' s <- 100
#' q <- 7
#' m <- s*2^q
#' obj <- SansSouciDyadic(m, s, flavor = "window.size", direction = "top-down")
#' 
#' mu <- gen.mu.leaves(m = m, K1 = 8, d = 0.9, grouped = TRUE, 
#'   setting = "const", barmu = 3, leaf_list = obj$input$leaves)
#' pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
#' 
#' alpha <- 0.05
#' res <- fit(obj, alpha, pvalues, "DKWM")
#' bound(res)
fit.SansSouciStruct <- function(object, alpha, p.values, 
                                family = c("DKWM", "HB", "trivial", "Simes", "Oracle"), 
                                flavor = c("tree", "partition"), ...) {
    family <- match.arg(family)
    flavor <- match.arg(flavor)
    if (family == "Oracle") {
        truth <- object$input$truth
        if (is.null(truth)) {
            stop("'truth' should be available for 'Oracle'. See ?SansSouci")
        }
    }
    params <- list(alpha = alpha, family = family, flavor = flavor)
    if (family %in% c("DKWM", "HB", "trivial")) {
        zeta <- switch(family,
                       "DKWM" = zeta.DKWM,
                       "HB" = zeta.HB,
                       "trivial" = zeta.trivial)
        struct <- object$input$struct
        leaves <- object$input$leaves
        if (flavor == "partition") {
            struct <- struct[length(struct)]
        }
        ZL <- zetas.tree(C = struct, leaf_list = leaves, 
                         method = zeta, pvalues = p.values, alpha = alpha)
        output <- list(p.values = p.values, ZL = ZL)
    } else if (family == "Simes") {
        thr <- SimesThresholdFamily(m, kMax =  m)(alpha) ## NB: we force K=m
        output <- list(p.values = p.values, thr = thr)
        class(object) <- "SansSouci" ## so that bound.SansSouci is called...
    } else if (family == "Oracle") {
        thr <- object$input$truth
        output <- list(p.values = p.values, thr = thr)
        class(object) <- "SansSouci" ## so that bound.SansSouci is called...
    }
    object$output <- output
    object$parameters <- params
    object
}

#' `label` Get the label of a post hoc method
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
label.SansSouciStruct <- function(object) {
    param <- object$parameters
    if (is.null(param)) {
        return(NULL)
    }
    lab <- param$family
    flv <- param$flavor
    
    if (!(lab %in% c("Simes", "Oracle"))) {
        lab <- sprintf("%s(%s)", lab, flv)
    } 
    return(lab)
}
pValues.SansSouciStruct <- function(object) {
    object$output$p.values
}

#' @export
bound.SansSouciStruct <- function(object, S = seq_len(nHyp(object)), 
                            what = c("TP", "FDP"), ...) {
    p.values <- pValues(object)
    if (max(S) > nHyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    struct <- object$input$struct
    leaves <- object$input$leaves
    ZL <- object$output$ZL
    
    if (object$param$flavor == "partition") {
        struct <- struct[length(struct)]
    }
    
    ord <- order(p.values)
    max_FP <- V.star(ord[S], C = struct, ZL, leaf_list = leaves)
    max_FP    
}
