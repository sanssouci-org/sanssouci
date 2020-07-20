#' Define a reference family
#' 
#' Define a reference family based on its name and parameter value
#'
#' @param refFamily A character value, the name of the reference family 
#' @param param A numeric value, the parameter value
#' @param digits An integer value, the number of decimal places to round to (by
#'   lower values)
#' 
#' @return A character value of the form "<refFamily>(<param>)"
#'   and 'value' is the numeric value of the
#'   parameter for this family, e.g. "Simes(0.05)". Currently, 'FamilyName'
#'   should be either "Simes" (or equivalenlty, "Linear") or "Beta".
#'
#' @seealso fromFamily
#' @aliases family
#' @export
#' @examples
#' 
#' alpha <- 0.1
#' fam <- list("Simes" = toFamily("Simes", alpha),
#'              "Simes + calibration" = toFamily("Simes", cal$lambda))
toFamily <- function(refFamily, param, digits = 2) {
    prec <- 10^digits
    thr <- floor(param*prec)/prec
    patt <- paste("%s(%.", digits + 2, "s)", sep = "")  ## ie "%s(%.4s)" for digits=2
    sprintf(patt, refFamily, thr)
}

#' Retrieve information from a reference family
#' 
#' Retrieve the name and parameter of a reference family
#' 
#' @details this function is the reciprocal of \code{\link{toFamily}}
#' @seealso toFamily
#' @aliases family
fromFamily <- function(fam) {
    patt <- "(.*)\\((.*)\\)"
    list(refFamily = gsub(patt, "\\1", fam),
         param = as.numeric(gsub(patt, "\\2", fam)))
}