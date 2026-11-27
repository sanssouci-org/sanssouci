#' Fit Joint Error Rate with k-Family Wise Error Rate
#'
#' @inheritParams fit.sanssouci_data object alpha adapt stepdown
#' @param method A string in \{"simes", "lambdacalibration", "heterosimes"\}.
#'                "simes" correponds to naive Simes template.
#'                "lambdacalibration" corresponds to lambda calibration with
#'                different template (see `template` argument)
#'                "heterosimes" : not currently implemented but in process.
#' @param n_permutatitons An integer value, the number of permutations to be performed.
#'  Defaults to 1000
#' @param K An integer value in `[1,m]`, the number of elements in the
#'   template. Defaults to m, the number of hypothesis.
#' @param template A character value, the name of a threshold family. Should be
#'   one of "Simes", "Beta", or "Oracle".
#'
#'   - Simes: The classical family of thresholds introduced by
#'   Simes (1986).
#'   This family yields JER control if the test statistics are positively
#'   dependent (PRDS) under H0.
#'
#'   - Beta: A family of thresholds that achieves marginal kFWER control under
#'   independence
#'
#'   - Oracle A family such that the associated bounds correspond to the true
#'   numbers/proportions of true/false positives. "truth" must be available in
#'   object$input$truth.
#'
#' @param p_cdfs argument for method = "heterosimes"
#'
#' @param max_steps_down A numeric value, the maximum number of steps down to
#'   perform. Defaults to 10 (but the algorithm generally converges in 1 or 2
#'   steps).
#' @param verbose A boolean, well you know verbose argument,
#' by default FALSE.
#'
#' @returns a ref_family_kfwer object.
#' @export
#'
#' @examples
#' # todo
fit_kfwer.sanssouci_data <- function(object, alpha,
                                     method = c(
                                       "simes",
                                       "lambdacalibration",
                                       "heterosimes"
                                     ),
                                     adapt = FALSE, stepdown = FALSE,
                                     n_permutations = 1000,
                                     K = n_hyp.SansSouci(object),
                                     template = c("Simes", "Beta", "Oracle",
                                                  "notip"),
                                     p_cdfs = NULL, max_stepdown = 10L,
                                     verbose = FALSE) {
  match.arg(method)
  match.arg(template)

  if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' should be a numeric value in [0,1]")
  }
  if (!is.null(p_cdfs) & method != "heterosimes") {
    message(
      "'p_cdfs' is not used in the case method != 'heterosimes'.",
      "Please use method = 'heterosimes'."
    )
    p_cdfs <- NULL
  }

  if (!stepdown) {
    max_stepdown <- 0L
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

    if (!is.null(template) && template != "Simes") {
      message("With method = 'simes', template is set to 'Simes'")
      template <- "Simes"
    }


    if (K != m) {
      message("For the method = 'simes', K is set to m.")
      K <- m
    }

    if (!adapt && !stepdown) {
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

    row_test_fun <- get_test(object) ## todo add get_test.sanssouci_data
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
    Y <- object$Y
    alternative <- object$alternative
    type <- object$type
    n_obs <- n_observations(object)


    ## core function
    if (!is.matrix(p_values)) {
      p_values <- matrix(p_values, nrow = 1)
      rownames(p_values) <- object$contrast_name # todo get_contrast_name ?
    }

    if (template %in% c("Simes", "Beta")) {
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
          categ = null_groups,
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
          max_steps_down = max_stepdown,
          piv_stat0 = NULL
        )

        thr <- calib$thr
        lambda <- calib$lambda

        if (verbose) {
          dt <- Sys.time() - t0
          cat("done (", format(dt), ")\n", sep = "")
        }
      } else if (n_permutations == 0) {
        if (template == "Beta" && K == m) {
          message("For the template = 'Beta',  we recommend choosing K < m.")
        }

        if (!adapt && !stepdown) {
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
        lambda <- NULL # ou ici alpha ?
      }
    } else if (template == "Oracle"){
      hat_m = m
      lambda = NULL
      thr = object$truth
    }
  } else if (method == "heterosimes") {
    stop(
      "For now, the heterosimes is not currently implemented.",
      "Come back soon ;)"
    )
  }

  ## create object
  object <- new_ref_family_kfwer(
    alpha = alpha, method = method,
    adapt = adapt, stepdown = stepdown,
    n_permutations = n_permutations, K = K,
    template = template, p_cdfs = p_cdfs,
    p_values = p_values, thr = thr,
    hat_m0 = hat_m, lambda = lambda
  )
}

new_ref_family_kfwer <- function(alpha, method, adapt, stepdown,
                                 n_permutations, K, template, p_cdfs,
                                 p_values, thr, hat_m0, lambda) {
  parameters <- list(
    alpha = alpha,
    setting = "kfwer",
    method = method,
    adapt = adapt,
    stepdown = stepdown,
    n_permutations = n_permutations,
    K = K,
    template = template,
    p_cdfs = p_cdfs
  )

  object <- list(
    parameters <- parameters,
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
