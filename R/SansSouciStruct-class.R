#' Create an object of class 'SansSouciStruct'
#' 
#' 
#' @param structure A (dyadic) structure
#' @param leaves A list of leaves
#' @param truth An optional numeric vector of $m$ values in ${0,1}$, the status of each null hypothesis (0 means H0 is true, 1 means H1 is true). Typically used in simulations.
#' @export
#' @examples
#' s <- 100
#' q <- 7
#' m <- s*2^q
#' 
#' dd <- dyadic.from.window.size(m, s, method = 2)
#' leaf_list <- dd$leaf_list
#' C <- dd$C
#' obj1 <- SansSouciStruct(C, leaf_list)
#' obj <- SansSouciDyadic(m, s, flavor = "window.size", direction = "top-down")
SansSouciStruct <- function(structure, leaves, truth = NULL) {
    input = list(structure = structure, 
                 leaves = leaves)
    input$truth <- truth
    obj <- structure(list(input = input,
                          parameters = NULL,
                          alpha = NULL,
                          output = NULL), 
                     class = "SansSouciStruct")
    obj
}

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
    SansSouciStruct(structure = dd$C, leaves = dd$leaf_list, ...)
}


#' `nHyp`: get the number of hypotheses
#' 
#' @rdname SansSouciStruct-methods
#' @param object An object of class `SansSouciStruct`
#' @export
nHyp.SansSouciStruct <- function(object) {
    length(unlist(object$input$leaves))
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
#' K1 <- 8
#' r <- 0.9
#' m1 <-  r*K1*s
#' barmu <- 3
#' mu <- gen.mu.leaves(m = m, K1 = K1, d = r, grouped = TRUE, 
#'   setting = "const", barmu = barmu, leaf_list = leaf_list)
#' pvalues <- gen.p.values(m = m, mu = mu, rho = 0)
#' alpha <- 0.05
#' res <- fit(obj, alpha, pvalues, "DKWM")
#' bound(res)
fit.SansSouciStruct <- function(object, alpha, p.values, 
                                zeta = c("DKWM", "HB", "trivial"), ...) {
    zeta <- match.arg(zeta)
    zetaFUN <- switch(zeta,
                      "DKWM" = zeta.DKWM,
                      "HB" = zeta.HB,
                      "trivial" = zeta.trivial)
    structure <- object$input$structure
    leaves <- object$input$leaves
    ZL <- zetas.tree(C = structure, leaf_list = leaves, 
                     method = zetaFUN, pvalues = p.values, alpha = alpha)
    object$output <- list(p.values = p.values, ZL = ZL)
    
    object$parameters <- list(
        alpha = alpha,
        zeta = zeta)
    
    object
}

label.SansSouciStruct <- function(object) {
    object$parameters$zeta
}

pValues.SansSouciStruct <- function(object) {
    object$output$p.values
}

#' @export
bound.SansSouciStruct <- function(object, S = seq_len(nLeaves(object)), 
                            what = c("TP", "FDP"), all = FALSE, ...) {
    p.values <- pValues(object)
    if (max(S) > nHyp(object)) {
        stop("'S' is not a subset of hypotheses")
    }
    lab <- label(object)

    structure <- object$input$structure
    leaves <- object$input$leaves
    
    ord <- order(p.values)
    ZL <- object$output$ZL
    
    max_FP <- maxFP_struct(p.values, S, ZL, structure, leaves, all = TRUE)
    bounds <- formatBounds(max_FP, idxs = seq_along(S), lab = lab, what = what, all = all)
    bounds
    
    return(bounds)
}


maxFP_struct <- function(p.values, S, ZL, C, leaf_list, 
                            what = c("TP", "FDP"), all = FALSE, ...) {
    if (max(S) > length(p.values)) {
        stop("'S' is not a subset of hypotheses")
    }
    ord <- order(p.values)

    V <- NULL
    if (!all) {
        V <- V.star(ord[S], C, ZL, leaf_list)
    } else {
        V <- sapply(S, FUN = function(ii) {
            V.star(ord[1:ii], C, ZL, leaf_list)
        })
    }
    return(V)
}

## TODO: implement method Oracle, Simes, + flavor tree/partition... and hybrid!