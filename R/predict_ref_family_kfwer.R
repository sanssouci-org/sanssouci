#' Post hoc confidence bounds on the true/false positives
#'
#' @param object An object of class 'sanssouci_data'
#' @param selection A subset of indices
#' @param what A character vector, the names of the post hoc bounds to be
#'   computed, among:
#'
#' - FP: Upper bound on the number of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the number of true positives in the 'x' most
#' significant items
#' - FDP: Upper bound on the proportion of false positives in the 'x' most
#' significant items
#' - TP: Lower bound on the proportion of true positives in the 'x' most
#' significant items
#'
#' Defaults to `c("TP", "FDP")`
#' @param all A logical value: should the bounds for all ordered subsets of `S`
#' be returned? If `FALSE` (the default), only the bound for `S` is returned
#' @param contrast_name A character value, the selected contrast.
#' Should be chosen in \code{x$input$contrast_name}.
#' @param ... Not used
#'
#' @return If `all` is `FALSE` (the default), only the value of the bound is
#' returned. Otherwise, a `data.frame` is return, with |S| rows and 4 columns:
#' - x: Number of most significant items selected
#' - label: Label for the procedure, typically of the form 'family(param)'
#' - bound: Value of the post hoc bound
#' - stat: Type of post hoc bound, as specified by argument `bound`.
#'
#' @importFrom stats predict
#' @importFrom rlang missing_arg is_missing
#' @export
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' table(groups)
#' a <- new_sanssouci_data(Y = expr_ALL, X = groups)
#'
#' res <- fit_kfwer(a, alpha = 0.1, method = "simes")
#' predict(res)
#'
#' # idem for all possible subsets (sorted by p-value)
#' bounds <- predict(res, all = TRUE)
#' head(bounds)
#'
#' # post hoc bound on a subset
#' S <- which(get_p_values(res) < 0.01)
#' predict(res, S)
predict.ref_family_kfwer <- function(object,
                              S = seq_len(get_n_hyp(object)),
                              what = c("TP", "FDP"), all = FALSE,
                              contrast_name = object$input$contrast_name, ...) {
  # S should be size of D
  if (!all(contrast_name %in% object$input$contrast_name)) {
    stop("contrast_name must be in object$input$contrast_name")
  }
  p_values <- get_p_values(object)
  thr <- get_thresholds(object)
  lab <- get_templates(object)

  # If sansSouci object is two sample test (1 contrast) /struct (0 contrast)
  if (length(object$input$contrast_name) <= 1) {
    bounds <- posthoc_bound(p_values,
                            S = S, thr = thr, lab = lab,
                            what = what, all = all
    )
    if (!all) {
      bounds <- bounds[, "bound"]
      if (length(bounds) > 1) {
        names(bounds) <- what
      }
    }

    return(bounds)
    # If there are more than 1 tested contrast
  } else {
    bounds_list <- list()
    for (contrasts in contrast_name) {
      # Only keep p-values of contrasts
      p_values_contrast <- p_values[contrasts, ]
      bounds <- posthoc_bound(p_values_contrast,
                              S = S, thr = thr, lab = lab,
                              what = what, all = all
      )
      if (!all) {
        bounds <- bounds[, "bound"]
        if (length(bounds) > 1) {
          names(bounds) <- what
        }
      }
      bounds_list[[contrasts]] <- bounds
    }
    # when only one contrast is display, a data.frame is returned
    if (length(contrast_name) == 1) {
      return(bounds_list[[1]])
    }
    return(bounds_list)
  }
}
