#' Create an object of class 'ref_family_levelsets'
#'
#' @param struct A list...
#' @param zetas_list A list...
#' @param alpha A value...
#' @param method A character value...
#' @param step_down A boolean value...
#' @param ... Other parameters to be listed.
#' @return An object of class `ref_family_levelsets`
#' @export
#'
new_ref_family_levelsets <- function(
    struct,
    zetas_list,
    alpha,
    method = c("Simes", "DKWM", "Bretagnolle", "KR", "Wellner"),
    step_down = FALSE,
    ...
) {
  method <- match.arg(method)
}

#' Fit Joint Error Rate with ...
#'
#' @export
fit_levelsets <- function(object, alpha, method, ...) UseMethod("fit_levelsets")


#' Fit Joint Error Rate with ...
#'
#' @inheritParams fit.sanssouci_data
#' @export

fit_levelsets.sanssouci_data <- function(
    object,
    alpha,
    method = c("Simes", "DKWM", "Bretagnolle", "KR", "Wellner"),
    adapt = FALSE,
    step_down = FALSE,
    ...
) {
  method <- match.arg(method)
  stop("Setting 'levelsets' not currently available")
}
  