fit_kwer.sanssouci_data <- function(object, alpha,
                                    method = c(
                                      "simes",
                                      "lambdacalibration",
                                      "heterosimes"
                                    ),
                                    adapt = FALSE, stepdown = FALSE,
                                    n_permutations = 1000,
                                    K = n_hyp.SansSouci(object),
                                    template = c("Simes", "Beta", "notip"),
                                    p_cdfs = NULL) {
  match.arg(method)
  match.arg(template)

  if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' should be a numeric value in [0,1]")
  }

  p_values = get_p_values(object)

  m = n_hyp.SansSouci(object)

  if (method == "simes") {
    if (!is.null(p_cdfs)) {
      message("'p_cdfs' is not used in the case.",
              "Please use method = 'heterosimes'.")
      p_cdfs = NULL
    }

    if (!is.null(n_permutations) & n_permutations != 0) {
      warning(
        "'n_permutations' is not used for method = 'simes'.",
        "To performe calibration, please use",
        "method = 'lambdacalibration'"
      )
      n_permutations = NULL
    }

    if(!is.null(template) & template != "Simes"){
      message("With method = 'simes', template is set to 'Simes'")
      template = "Simes"
    }


    if(K != m){
      message("For the Simes method, K is set to m.")
      K = m
    }

    if(!adapt & !stepdown){
      hat_m = m
    }

    thr = simes$t(alpha, 1:hat_m, hat_m)
    lambda = NULL #ou ici alpha ?

  } else if (method == "lambdacalibration") {

    ## check parameters
    if (template == "notip") {
      stop("'notip' template is not currently implemented. Come back soon ;)")
    }





    ## core function
    if (!is.matrix(cal$p_value)) {
      cal$p_value <- matrix(cal$p_value, nrow = L)
      rownames(cal$p_value) <- object$input$contrast_name
    }
    if (!is.null(cal$estimate) && !is.matrix(cal$estimate)) {
      cal$estimate <- matrix(cal$estimate, nrow = L)
      rownames(cal$estimate) <- object$input$contrast_name
    }
    p_values <- cal$p_value
    if (B > 0 && template %in% c("Simes", "Beta")) {
      t0 <- Sys.time()
      if (verbose) {
        cat("Randomization p-values...")
      }
      if (do_p0) {
        null_groups <- NULL
        if (type == "1 sample") { # sign-flipping
          null_groups <- replicate(B, rbinom(n, 1, 0.5) * 2 - 1)
        } else if (type == "2 samples") { # permutation
          null_groups <- replicate(B, sample(groups))
        } else if (type == "Continuous case, correlation test") {
          # continuous covariate
          null_groups <- replicate(B, sample(groups))
        } else if (type == "linear model") {
          null_groups <- replicate(B, sample(1:n, n, replace = TRUE))
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
      } else if (verbose) {
        cat("skipped computation (already done)\n")
      }
      t0 <- Sys.time()
      if (verbose) {
        cat("Calibration...")
      }
      calib <- sanssouci::calibrate(
        p0 = p0, m = m, alpha = alpha,
        family = family, K = K,
        p = p_values,
        max_steps_down = max_steps_down,
        piv_stat0 = piv_stat0
      )
      if (verbose) {
        dt <- Sys.time() - t0
        cat("done (", format(dt), ")\n", sep = "")
      }
    }

  } else if (method == "heterosimes") {
    stop(
      "For now, the heterosimes is not currently implemented.",
      "Come back soon ;)"
    )
  }

  ## create object
  object <- new_ref_family_kfwer(alpha = alpha, method = method,
                                 adapt = adapt, stepdown = stepdown,
                                 n_permutations = n_permutations, K = K,
                                 template = template, p_cdfs = p_cdfs,
                                 p_values = p_values, thr = thr,
                                 hat_m0 = hat_m0, lambda = lambda)


}

new_ref_family_kfwer <- function(alpha, method, adapt, stepdown,
                                   n_permutations, K, template, p_cdfs,
                                   p_values, thr, hat_m0, lambda){

  parameters <- list(
    alpha = alpha,
    setting = "kfwer",
    method = method,
    adapt = adapt,
    stepdown = stepdown,
    n_permutations = n_permutations,
    K = K,
    template = templates,
    p_cdfs = p_cdfs
  )

  object <- list(
    parameters <- parameters,
    p_values = p_values,
    thr = thr,
    hat_m0 = hat_m0,
    lambda = lambda
  )

  object <- structured(
    object,
    class = "ref_family_kfwer"
  )

}
