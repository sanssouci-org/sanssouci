#' Fit Joint Error Rate with k-Family Wise Error Rate
#'
#' @export
fit_kfwer <- function(object, alpha, method, ...) UseMethod("fit_kfwer")

#' Fit Joint Error Rate with k-Family Wise Error Rate
#'
#' @inheritParams fit.sanssouci_data object alpha adapt step_down
#' @param method A string in \{"Simes", "lambda-calibration", "Het Simes"\}.
#'   "Simes" corresponds to the Simes reference family. "lambda-calibration"
#'   corresponds to lambda calibration with different template (see `template`
#'   argument) "Het Simes": not currently implemented.
#' @param n_permutatitons An integer value, the number of permutations to be
#'   performed. Defaults to 1000.
#' @param K An integer value in `[1,m]`, the number of elements in the template.
#'   Defaults to m, the number of hypothesis.
#' @param template A character value, the name of a threshold family. Should be
#'   one of "Simes", "Beta", or "Oracle".
#'
#'   - Simes: The classical family of thresholds introduced by
#'   Simes (1986). This family yields JER control if the test statistics are
#'   positively dependent (PRDS) under H0.
#'
#'   - Beta: A family of thresholds that achieves marginal kFWER control under
#'   independence
#'
#'   - Oracle A family such that the associated bounds correspond to the true
#'   numbers/proportions of true/false positives. "truth" must be available in
#'   object$truth.
#'
#' @param p_cdfs argument for method = "Heterogeneous Simes"
#'
#' @param max_steps_down A numeric value, the maximum number of steps down to
#'   perform. Defaults to 10 (but the algorithm generally converges in 1 or 2
#'   steps).
#' @param verbose A boolean, well you know verbose argument, by default FALSE.
#'
#' @returns a ref_family_kfwer object.
#' @export
#'
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' table(groups)
#' a <- new_sanssouci_data(Y = expr_ALL, X = groups)
#'
#' res <- fit_kfwer(a, alpha = 0.1, method = "Simes")
#' res <- fit_kfwer(a, alpha = 0.1, method = "lambda-calibration", step_down = TRUE)
fit_kfwer.sanssouci_data <- function(object, alpha,
                                     Method = c(
                                       "Simes",
                                       "lambda-calibration",
                                       "Heterogeneous Simes"
                                     ),
                                     adapt = FALSE, step_down = FALSE,
                                     n_permutations = NULL,
                                     K = get_n_hyp(object),
                                     template = c("Simes", "Beta", "Oracle",
                                                  "Notip"),
                                     p_cdfs = NULL, max_step_down = 10L,
                                     verbose = FALSE) {
  method <- match.arg(method)
  template <- match.arg(template)

  print(method)
  print(template)

  if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' should be a numeric value in [0,1]")
  }
  if (!is.null(p_cdfs) & method != "Heterogeneous Simes") {
    message("'p_cdfs' only used for method 'Heterogeneous Simes'.")
    p_cdfs <- NULL
  }

  if (!step_down) {
    max_step_down <- 0L
  }

  p_values <- get_p_values(object)

  m <- get_n_hyp(object)



  if (method == "simes") {

    if (!is.null(n_permutations) && n_permutations != 0) {
      warning(
        "'n_permutations' is not used for method = 'simes'.",
        "To performe calibration, please use",
        "method = 'lambdacalibration'"
      )
      n_permutations <- NULL
    }

    if ( template != "Simes") {
      message("With method = 'simes', template is set to 'Simes'")
      template <- "Simes"
    }


    if (K != m) {
      message("For the method = 'simes', K is set to m.")
      K <- m
    }

    if (!adapt && !step_down) {
      hat_m <- m
    } else {
      ## todo
      hat_m <- m
    }

    thr <- simes$t(alpha, 1:hat_m, hat_m)
    lambda <- NULL # ou ici alpha ?
  } else if (method == "lambdacalibration") {
    ## check parameters
    if (template == "notip") {
      stop("'notip' template is not currently implemented. Come back soon ;)")
    }

    row_test_fun <- get_test_fun(object)
    if (is.null(row_test_fun)) {
      stop(
        "Test function must be refered in new_sanssouci_data",
        "to compute lambda calibration.",
        "If you want method = 'lambdacalibration', please refer 'test' in",
        "new_sanssouci_data.",
        "If you want to use your own p-values, please use method = 'simes'",
        "in fit_kfwer."
      )
    }

    ## getter on object
    Y <- get_y(object)
    alternative <- object$alternative
    type <- get_test_type(object)
    n_obs <- get_n_observations(object)
    groups <- get_groups_to_permute(object)


    ## core function
    if (!is.matrix(p_values)) {
      p_values <- matrix(p_values, nrow = 1)
      rownames(p_values) <- object$contrast_name # todo get_contrast_name ?
    }

    if (template %in% c("Simes", "Beta")) {
      if(is.null(n_permutations)){
        n_permutations <- 1000
        message("for lambadacalibration, if _n_permutations is null, then ",
                "n_permutations is set to 1000 by default")
      }
      if (n_permutations > 0) {
        t0 <- Sys.time()
        if (verbose) {
          cat("Randomization p-values...")
        }
        null_groups <- NULL
        if (type == "1 sample") { # sign-flipping
          null_groups <- replicate(n_permutations,
                                   rbinom(n_obs, 1, 0.5) * 2 - 1)
        } else if (type == "2 samples") { # permutation
          null_groups <- replicate(n_permutations, sample(groups)) ## todo
        } else if (type == "Continuous case, correlation test") {
          # continuous covariate
          null_groups <- replicate(n_permutations, sample(groups)) ## todo
        } else if (type == "linear model") {
          null_groups <- replicate(n_permutations,
                                   sample(1:n_obs, n_obs, replace = TRUE))
        } else {
          stop("type is not defined please put an issue.")
        }
        p0 <- minibatch_row_test_fun(
          row_test_fun = row_test_fun, Y = Y,
          groups = null_groups,
          alternative = alternative,
          max_batch_size = 1e6, m = m
        )

        if (verbose) {
          dt <- Sys.time() - t0
          cat("done (", format(dt), ")\n", sep = "")
        }
        t0 <- Sys.time()
        if (verbose) {
          cat("Calibration...")
        }
        calib <- calibrate(
          p0 = p0, m = m, alpha = alpha,
          family = template, K = K,
          p = p_values,
          max_steps_down = max_step_down,
          piv_stat0 = NULL
        )

        thr <- calib$thr
        lambda <- calib$lambda
        hat_m <- m

        if (verbose) {
          dt <- Sys.time() - t0
          cat("done (", format(dt), ")\n", sep = "")
        }
      } else if (n_permutations == 0) {
        if (template == "Beta" && K == m) {
          message("For the template = 'Beta',  we recommend choosing K < m.")
        }

        if (!adapt && !step_down) {
          hat_m <- m
        } else {
          ## todo
          hat_m <- m
        }

        obj_template <- switch(template,
          "Simes" = simes,
          "Beta" = beta
        )

        thr <- obj_template$t(alpha, 1:hat_m, hat_m)
        lambda <- alpha
      }
    } else if (template == "Oracle"){
      
      if(is.null(get_truth(object))){
        stop("Please provide truth in new_sanssouci_data.")
      }
      hat_m <- m
      lambda <- NULL
      thr <- get_truth(object)
    }
  } else if (method == "Heterogeneous Simes") {
    stop("'Heterogeneous Simes' method not implemented yet")
  }

  ## create object
  object <- new_ref_family_kfwer(
    alpha = alpha, 
    method = method,
    adapt = adapt, 
    step_down = step_down,
    n_permutations = n_permutations, 
    K = K,
    template = template, 
    p_cdfs = p_cdfs,
    p_values = p_values, 
    thr = thr,
    hat_m0 = hat_m, 
    lambda = lambda
  )
  return(object)
}



new_ref_family_kfwer <- function(alpha, method, adapt, step_down,
                                 n_permutations, K, template, p_cdfs,
                                 p_values, thr, hat_m0, lambda) {
  parameters <- list(
    setting = "kFWER",
    alpha = alpha,
    method = method,
    adapt = adapt,
    step_down = step_down,
    n_permutations = n_permutations,
    K = K,
    template = template,
    p_cdfs = p_cdfs
  )

  object <- list(
    parameters = parameters,
    p_values = p_values,
    thr = thr,
    hat_m0 = hat_m0,
    lambda = lambda
  )

  object <- structure(
    object,
    class = "ref_family_kfwer"
  )
}


#' `get_p_values`: get p-values
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_p_values.ref_family_kfwer <- function(object) {
  object$p_values
}

get_thresholds <- function(object) UseMethod("get_thresholds")
#' `get_thresholds`: get p-values
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_thresholds.ref_family_kfwer <- function(object) {
  object$p_values
}

get_templates <- function(object) UseMethod("get_templates")
#' `get_templates`: get p-values
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_templates.ref_family_kfwer <- function(object) {
  object$parameters$template
}
