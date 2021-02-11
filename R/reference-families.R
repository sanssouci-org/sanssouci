t_inv_linear <- function(y, k, m) {
    y * m / k;
}

t_linear <- function(alpha, k, m) {
    alpha * k / m;
}


t_inv_beta <- function(y, k, m) {
    pbeta(y, k, m + 1 - k);
}

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