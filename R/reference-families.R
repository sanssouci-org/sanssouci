#' Reference families
#'
#' @param alpha numeric value in (0,1)
#' @param y numeric value
#' @param k integer value in `[1,m]`
#' @param m integer value
#' @name reference-families
#' @examples
#' m <- 10
#' alpha <- 0.05
#' thr <- t_linear(alpha, 1:m, m)
#' t_inv_linear(thr[3], 3, m)
#' 
#' thr <- t_beta(alpha, 1:m, m)
#' t_inv_beta(thr[3], 3, m)
NULL
#> NULL

#' @name reference-families
#' @export
t_inv_linear <- function(y, k, m) {
    y * m / k;
}

#' @name reference-families
#' @export
t_linear <- function(alpha, k, m) {
    alpha * k / m;
}

#' @name reference-families
#' @importFrom stats pbeta
#' @export
t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}

#' @name reference-families
#' @importFrom stats qbeta
#' @export
t_beta <- function(alpha, k, m) {
    qbeta(alpha, k, m + 1 - k);
}

Simes <- list(name = "Simes",
              t = function(lambda, k, m) lambda * k / m,
              t_inv = function(y, k, m) y * m / k)

Beta <- list(name = "Beta", 
            t = function(lambda, k, m) qbeta(lambda, 1:m, m + 1 - k),
            t_inv = function(y, k, m) pbeta(y, k, m + 1 - k))

check_ref_fam <- function(fam) {
    ## check types
    stopifnot(typeof(fam) == "list")
    
    ## check names
    nms <- c("name", "t", "t_inv") 
    for (nm in nms) {
        if (is.null(fam[[nm]])) {
            stop("Element '", nm, "' missing")
        }
    }
    
    ## check types
    name <- fam[["name"]]
    if (is.null(name)) {
        stop("Element 'name' missing")
    } else {
        stopifnot(typeof(name) == "character")
    }
    
    FUN <- fam[["t"]]
    if (is.null(FUN)) {
        stop("Element 't' missing")
    } else {
        stopifnot(mode(FUN) == "function")
    }

    FUN <- fam[["t_inv"]]
    if (is.null(FUN)) {
        stop("Element 't_inv' missing")
    } else {
        stopifnot(mode(FUN) == "function")
    }
}