#' Create an object of class 'SansSouci'
#'
#' @param Y A matrix of \eqn{D} variables (hypotheses) by \eqn{n} observations
#' @param X A design matrix of \eqn{n} observations by \eqn{p} explanatory
#'   variables in the linear model
#' @param contrast A contrast matrix of size \eqn{L} tested contrasts (in row)
#'   and \eqn{p} columns corresponding to the parameters to be tested
#' @param groups A numeric vector of \eqn{n} values in \eqn{0, 1}, the groups of
#'   observations on which to perform two-sample tests
#' @param truth An optional numeric vector of \eqn{m} values in \eqn{0, 1},
#'   the status of each null hypothesis (0 means H0 is true, 1 means H1 is
#'   true). Typically used in simulations.
#' @return An object of class "SansSouci"
#' @export
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' table(groups)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' plot(res, xmax = 500)
#' res_beta <- fit(res, B = 100, alpha = 0.1, family = "Beta", K = 10)
#' plot(res_beta, xmax = 500)
#'
SansSouci <- function(Y, X = NULL, contrast = NULL, groups = NULL,
                      truth = NULL) {
  # Are the arguments present?
  signature <- c(!is.null(X), !is.null(contrast), !is.null(groups))


  n <- ncol(Y)
  D <- nrow(Y)

  # None arguments: one sample case
  if (identical(signature, c(FALSE, FALSE, FALSE))) {
    groups <- rep(1, ncol(Y))
    type <- "1 sample"
    L <- 1
    P <- 1
    contrast_name <- "one sample test"
    ## Only groups is given: detection of one, two samples or continuous
  } else if (identical(signature, c(FALSE, FALSE, TRUE))) {
    ugroups <- unique(groups)
    n_groups <- length(ugroups)

    if (n_groups > 1) {
      categ_check(groups, ncol(Y))
    }
    if (n_groups == 1) {
      type <- "1 sample"
      contrast_name <- "One sample test"
    } else if (n_groups == 2) {
      type <- "2 samples"
      contrast_name <- "Group 1 vs group 0"
      # } else if (n_groups == n) {
    } else { # by default consider a continuous
      type <- "Continuous case, correlation test"
      contrast_name <- "Correlation test"
      if (n_groups != n) {
        message(
          "Ties were detected in the grouping variable, but a continuous",
          " vector was provided for groups, and a correlation test is being",
          " performed. If you intend to test associations across more than",
          " two groups, please provide a full design matrix using the X",
          " argument instead."
        )
      }
    }

    L <- 1
    P <- 1

    ## X and contrast are given: linear model case
  } else if (identical(signature, c(TRUE, TRUE, FALSE))) {
    .check_lm_test(t(Y), X, contrast)
    type <- "linear model"
    L <- nrow(contrast)
    P <- ncol(X)
    if (is.null(rownames(contrast))) {
      contrast_name <- paste("contrast", seq_len(L))
      rownames(contrast) <- contrast_name
    } else {
      contrast_name <- rownames(contrast)
    }
    ## The contrast matrix is not provided
  } else if (identical(signature, c(TRUE, FALSE, FALSE))) {
    stop("Please give a contrast matrix in `contrast`")
    ## The design matrix is not provided
  } else if (identical(signature, c(FALSE, TRUE, FALSE))) {
    stop("Please give a design matrix in `X`")
    ## The vector groups is provided with the contrast or the design matrix
  } else if (identical(signature, c(TRUE, FALSE, TRUE)) ||
               identical(signature, c(FALSE, TRUE, TRUE)) ||
               identical(signature, c(TRUE, TRUE, TRUE))) {
    stop("Please provide either an `X` design matrix and a `contrast`  matrix,
         or only a `groups` design vector.")
  }

  if (!is.null(truth)) {
    all_0 <- identical(truth, rep(0, nrow(Y)))
    all_1 <- identical(truth, rep(1, nrow(Y)))
    if (!all_0 && !all_1) {
      categ_check(truth, nrow(Y))
    }
  }
  input <- list(
    Y = Y,
    X = X,
    C = contrast,
    groups = groups,
    n_obs = n,
    n_contrasts = L,
    n_dimensions = D,
    n_variables = P,
    contrast_name = contrast_name,
    m = D * L,
    type = type
  )
  input$truth <- truth
  obj <- structure(
    list(
      input = input,
      parameters = NULL,
      output = NULL
    ),
    class = "SansSouci"
  )
  obj
}


#' Create an object of class 'SansSouci' from simulations
#'
#' Create a 'SansSouci' object from simulation in the Gaussian
#' equi-correlated model
#'
#' @param ... Parameters to be passed to [sim_gaussian_equi]
#' @seealso [sim_gaussian_equi]
#' @export
#' @examples
#' obj <- SansSouciSim(
#'   m = 543, rho = 0.4, n = 210,
#'   pi0 = 0.8, SNR = 3, prob = 0.5
#' )
#' alpha <- 0.1
#'
#' # Adaptive Simes (lambda-calibration)
#' set.seed(542)
#' res <- fit(obj, B = 100, alpha = alpha, family = "Simes")
#' plot(res)
#' volcano_plot(res, q = 0.05, r = 0.2)
#'
#' # upper bound on number of signals if the entire data set
#' # (and corresponding lower bound on FDP)
#' predict(res)
#'
#' # confidence curve
#' plot(res)
#'
#' # comparison to other confidence curves
#' # Parametric Simes (no calibration -- assume positive dependence (PRDS))
#' res0 <- fit(obj, B = 0, alpha = alpha, family = "Simes")
#' res0
#'
#' # Oracle
#' oracle <- fit(obj, alpha = alpha, family = "Oracle")
#' oracle
#'
#' confs <- list(
#'   Simes = predict(res0, all = TRUE),
#'   "Simes+calibration" = predict(res, all = TRUE),
#'   "Oracle" = predict(oracle, all = TRUE)
#' )
#' plot_conf_curve(confs)
#'
#' \dontrun{
#' # Use wilcoxon tests instead of Welch tests
#' res <- fit(obj, B = 100, alpha = 0.1, row_test_fun = row_wilcoxon_tests)
#' volcano_plot(res, q = 0.05, r = 0.2)
#' }
#'
SansSouciSim <- function(...) {
  sim <- sim_gaussian_equi(...)
  SansSouci(Y = sim$X, groups = sim$categ, truth = sim$H)
}

#' Basic methods for class `SansSouci`
#'
#' @param object An object of class \code{SansSouci}
#' @name SansSouci-methods
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' table(groups)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#' print(a)
#' n_hyp(a)
#' n_observations(a)
#' label(a)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' label(res)
#' print(res)
#' str(get_p_values(res))
#' str(fold_changes(res))
#' str(thresholds(res))
#' vp <- volcano_plot(res, q = 0.05, r = 0.5)
#'
NULL
#> NULL

#' Generic functions for S3 class SansSouci
#'
#' @name all-generics
#' @param object An object. See individual methods for specifics
NULL
#> NULL

#' @describeIn all-generics ' Get the number of hypotheses
#' @param object An object. See individual methods for specifics
#' @export
n_hyp <- function(object) UseMethod("n_hyp")

#' `n_hyp`: get the number of hypotheses
#'
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
n_hyp.SansSouci <- function(object) {
  object$input$m
}

#' @describeIn all-generics Get the number of observations
#' @param object An object. See individual methods for specifics
#' @export
n_observations <- function(object) UseMethod("n_observations")

#' `n_observations` Get the number of observations
#'
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
n_observations.SansSouci <- function(object) {
  object$input$n_obs
}

#' @describeIn all-generics Get the label of hypotheses tested
#' @param object An object. See individual methods for specifics
#' @export
label <- function(object) UseMethod("label")

#' `label` Get the label of a post hoc method
#'
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
label.SansSouci <- function(object) {
  param <- object$parameters
  if (is.null(param)) {
    return(NULL)
  }
  lab <- param$family
  if (!(lab %in% c("Simes", "Oracle")) && (!is.null(param$K))) {
    lab <- sprintf("%s(K=%d)", lab, param$K)
  }
  return(lab)
}

#' Print 'SansSouci' objects
#'
#' @rdname SansSouci-methods
#' @param x An object of class `SansSouci`
#' @param ... Not used
#' @param verbose Should detailed output be printed? Defaults to FALSE
#' @importFrom utils str
#' @export
print.SansSouci <- function(x, ..., verbose = FALSE) {
  object <- x
  rm(x)
  cat("'SansSouci' object:\n")
  input <- object$input
  if (!is.null(input)) {
    cat("\tNumber of hypotheses:\t", n_hyp(object), "\n")
    if (!is.null(n_observations(object))) {
      cat("\tNumber of observations:\t", n_observations(object), "\n")
    }
    cat("\t", input$type, "\n", sep = "")
    cat("\n")
    if (verbose) {
      cat("Data:")
      cat("\n")
      str(input$Y)
      cat("\n")
    }
    truth <- input$truth
    if (!is.null(truth)) {
      cat("Truth: ")
      cat("\n")
      cat("\t", sum(truth), "false null hypotheses (signals)")
      pi0 <- 1 - mean(truth)
      cat(" out of ", n_hyp(object), " (pi0=", round(pi0, 3), ")", sep = "")
      cat("\n")
    }
  }
  params <- object$parameters
  if (!is.null(params)) {
    cat("Parameters:", "\n")
    cat("\tTest function:\t\t", params$fun_name, "\n", sep = "")
    cat("\tNumber of permutations:\tB=", params$B, "\n", sep = "")
    cat("\tSignificance level:\talpha=", params$alpha, "\n", sep = "")
    cat("\tReference family:\t", params$fam, "\n", sep = "")
    cat("\t\t(of size:\tK=", params$K, ")", "\n", sep = "")
    cat("\n")
  }
  output <- object$output
  if (!is.null(output)) {
    cat("Output:\n")
    if (verbose) {
      cat("\tp-values:\n")
      print(summary(output$p_value))
      cat("\tCalibrated thresholds:")
      if (all(output$thr %in% c(0, 1))) {
        print(table(output$thr))
        cat("\n")
      } else {
        print(summary(output$thr))
        cat("\tPivotal statistic:\n")
        print(summary(output$piv_stat))
      }
      cat("\n")
    }
    cat("\tCalibration parameter:\tlambda=", output$lambda, "\n", sep = "")
  }
  invisible(object)
}

#' @importFrom generics fit
#' @export
generics::fit

#' Fit SansSouci object
#' @param object An object of class `SansSouci`
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#' @param B An integer value, the number of permutations to be performed.
#'  Defaults to 1000
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of "two.sided" (default), "greater" or "less".
#' @param row_test_fun A (vectorized) test function.
#'  Defaults to [row_welch_tests].
#' @param family A character value, the name of a threshold family. Should be
#'   one of "Linear", "Beta" and "Simes", or "Oracle". "Linear" and "Simes"
#'   families are identical.
#'
#'   - Simes/Linear: The classical family of thresholds introduced by
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
#' @param max_steps_down A numeric value, the maximum number of steps down to
#'   perform. Defaults to 10 (but the algorithm generally converges in 1 or 2
#'   steps).
#' @param K An integer value in `[1,m]`, the number of elements in the
#'   reference family. Defaults to m
#' @param force A boolean value: should the permutation p-values and pivotal
#' statistics be re-calculated ? Defaults to `FALSE`
#' @param verbose A boolean value: should extra info be printed? Defaults to
#' `FALSE`
#' @param ... Not used
#' @return A 'fitted' object of class 'SansSouci'. It is a list of three
#'    elements
#'  - input: see [SansSouci]
#'  - param: the input parameters, given as a list
#'  - output: the outputs of the calibration, see [calibrate]
#' @export
#' @examples
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(
#'   m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
#'   prob = 0.5
#' )
#' res <- fit(obj, B = 100, alpha = 0.1)
#'
#' # confidence curve
#' plot(res)
#'
#' # confidence curve for a subset
#' S <- which(get_p_values(res) < 0.1 & get_fold_changes(res) > 0.3)
#' plot(res, S = S)
#'
#' # plot two confidence curves
#' res_beta <- fit(res, B = 100, alpha = 0.1, family = "Beta", K = 20)
#'
#' resList <- list("Linear" = res, "Beta" = res_beta)
#' bounds <- lapply(resList, predict, all = TRUE)
#' plot_conf_curve(bounds, xmax = 200)
#'
fit.SansSouci <- function(object, alpha, B = 1e3,
                          row_test_fun = NULL,
                          alternative = c("two.sided", "less", "greater"),
                          family = c("Simes", "Linear", "Beta", "Oracle"),
                          max_steps_down = 10L, K = n_hyp(object),
                          force = FALSE,
                          verbose = FALSE, ...) {
  alternative <- match.arg(alternative)
  family <- match.arg(family)
  if (family == "Oracle") {
    alpha <- NA_real_ ## alpha is not used by Oracle
    truth <- object$input$truth
    if (is.null(truth)) {
      stop("'truth' should be available for 'Oracle'. See ?SansSouci")
    }
  } else {
    if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
      stop("Argument 'alpha' should be a numeric value in [0,1]")
    }
  }
  Y <- object$input$Y
  X <- object$input$X
  C <- object$input$C
  groups <- object$input$groups
  m <- n_hyp(object)
  n <- n_observations(object)
  L <- object$input$n_contrasts
  fun_name <- NA_character_
  type <- object$input$type

  if (is.null(row_test_fun)) {
    if (type == "1 sample") {
      row_test_fun <- row_z_tests
      fun_name <- "row_z_tests"
    } else if (type == "2 samples") {
      row_test_fun <- row_welch_tests
      fun_name <- "row_welch_tests"
    } else if (type == "Continuous case, correlation test") {
      row_test_fun <- row_pearson_correlation_tests
      fun_name <- "row_pearson_correlation_tests"
    } else if (type == "linear model") {
      row_lm_test_with_x_c <- function(...) {
        row_lm_test(..., X = X, C = C)
      }
      row_test_fun <- row_lm_test_with_x_c
      fun_name <- "row_lm_test"
      groups <- matrix(1:n, ncol = 1)
    }
  } else {
    fun_name <- as.character(substitute(row_test_fun))
  }

  ## should we re-calculate p0?
  params <- object$parameters
  p0 <- object$output$p0
  do_p0 <- TRUE
  if (!is.null(params) && !force) {
    if (!is.null(p0)) {
      cond_n_permutation <- (params$B == B)
      cond_test_function <- all.equal(params$row_test_fun, row_test_fun)
      cond_alternative <- (params$alternative == alternative)
      do_p0 <- (!cond_n_permutation) || (!cond_test_function) ||
        (!cond_alternative)
    }
  }

  ## should we re-calculate the (first) pivotal statistic ?
  params <- object$parameters
  piv_stat0 <- NULL
  if (!is.null(params) && !force) {
    cond_family <- (params$family == family)
    cond_k <- (params$K == K)
    cond_p <- (!do_p0)
    if (cond_family && cond_k && cond_p) {
      piv_stat0 <- object$output$piv_stat
    }
  }
  object$parameters <- list(
    alpha = alpha,
    B = B,
    alternative = alternative,
    row_test_fun = row_test_fun,
    fun_name = fun_name,
    type = type,
    family = family,
    max_steps_down = max_steps_down,
    K = K
  )

  if (family == "Beta" && K == m) {
    warning("For the 'Beta' family we recommend choosing K < m")
  }

  cal <- row_test_fun(Y, groups, alternative = alternative)
  if (!is.matrix(cal$p_value)) {
    cal$p_value <- matrix(cal$p_value, nrow = L)
    rownames(cal$p_value) <- object$input$contrast_name
  }
  if (!is.null(cal$estimate) && !is.matrix(cal$estimate)) {
    cal$estimate <- matrix(cal$estimate, nrow = L)
    rownames(cal$estimate) <- object$input$contrast_name
  }
  p_values <- cal$p_value
  if (B > 0 && family != "Oracle") {
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

    cal$p0 <- p0
    cal <- c(cal, calib)
  } else { # no calibration!
    cal$lambda <- alpha
    thr <- switch(family,
      "Linear" = t_linear(alpha, seq_len(m), m),
      "Simes" = t_linear(alpha, seq_len(m), m),
      "Beta" = t_beta(alpha, seq_len(m), m),
      "Oracle" = object$input$truth
    )
    cal$thr <- thr
  }
  # }
  object$output <- cal
  object
}

#' @describeIn all-generics Get p-values
#' @param object An object. See individual methods for specifics
#' @export
get_p_values <- function(object) UseMethod("get_p_values")

#' `get_p_values`: get p-values
#'
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
get_p_values.SansSouci <- function(object) {
  object$output$p_value
}

#' @describeIn all-generics Get fold changes
#' @param object An object. See individual methods for specifics
#' @export
get_fold_changes <- function(object) UseMethod("get_fold_changes")

#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
get_fold_changes.SansSouci <- function(object) {
  object$output$estimate
}

## An attempt with 'setters'
# 'fold_changes<-' <- function(object, value) UseMethod("fold_changes")
# 'fold_changes<-.SansSouci' <- function(object, value) {
#     object$output$estimate <- value
#     object
# }


#' @describeIn all-generics Get thresholds
#' @param object An object. See individual methods for specifics
#' @export
thresholds <- function(object) UseMethod("thresholds")

#' `thresholds`: get thresholds
#'
#' @rdname SansSouci-methods
#' @param object An object of class `SansSouci`
#' @export
thresholds.SansSouci <- function(object) object$output$thr

#' Plot confidence bound on the true/false positives among most significant
#' items
#'
#' @param x An object of class 'SansSouci'
#' @param y Not used
#' @param xmax Right limit of the plot
#' @param ... Further arguments to be passed to \code{bound}
#' @export
#' @examples
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(
#'   m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
#'   prob = 0.5
#' )
#' res <- fit(obj, B = 100, alpha = 0.1)
#'
#' # confidence curve
#' plot(res)
#'
#' # confidence curve for a subset of hypotheses
#' S <- which(get_p_values(res) < 0.1 & get_fold_changes(res) > 0.3)
#' plot(res, S = S)
plot.SansSouci <- function(x, y, xmax = n_hyp(x), ...) {
  cb <- predict(x, all = TRUE, ...)
  plot_conf_curve(cb, xmax = xmax)
}

#' Post hoc confidence bounds on the true/false positives
#'
#' @param object An object of class 'SansSouci'
#' @param S A subset of indices
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
#'
#' # Generate Gaussian data and perform multiple tests
#' obj <- SansSouciSim(
#'   m = 502, rho = 0.5, n = 100, pi0 = 0.8, SNR = 3,
#'   prob = 0.5
#' )
#' res <- fit(obj, B = 100, alpha = 0.1)
#'
#' # post hoc bound on the set of all hypotheses
#' predict(res)
#'
#' # idem for all possible subsets (sorted by p-value)
#' bounds <- predict(res, all = TRUE)
#' head(bounds)
#'
#' # post hoc bound on a subset
#' S <- which(get_p_values(res) < 0.01)
#' predict(res, S)
predict.SansSouci <- function(object,
                              S = seq_len(min(
                                n_hyp(object),
                                object$input$n_dimensions
                              )),
                              what = c("TP", "FDP"), all = FALSE,
                              contrast_name = object$input$contrast_name, ...) {
  # S should be size of D
  if (!all(contrast_name %in% object$input$contrast_name)) {
    stop("contrast_name must be in object$input$contrast_name")
  }
  p_values <- get_p_values(object)
  thr <- thresholds(object)
  lab <- label(object)

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

#' Volcano plot for a `SansSouci` object
#'
#' @param fold_changes An optional vector of fold changes, of the same length
#' as `n_hyp(x)`, used for volcano_plot x-axis.
#' If not specified, `get_fold_changes(x)`
#' is used.
#' @param p_values A vector of p-values, of the same length as `n_hyp(x)`,
#' used for volcano_plot y-axis. If not specified, `get_p_values(x)` is used
#' @param contrast_name A character value, the selected contrast. Should be
#' chosen in `x$input$contrast_name`.
#' @inheritParams volcano_plot.numeric
#'
#' @details The default is to use the fold changes and p-values from the input
#' SansSouci object. See the vignette
#'   <https://sanssouci-org.github.io/sanssouci/articles/post-hoc_differential-expression_RNAseq.html#custom-statistics-example-using-limma-voom> # nolint
#'   for an example where custom fold changes and p-values are used.
#'
#' @export
#'
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' volcano_plot(res, q = 0.2, r = 0.2, ylim = c(0, 4), feature_label = "gene")
volcano_plot.SansSouci <- function(x,
                                   contrast_name = x$input$contrast_name[1],
                                   fold_change = get_fold_changes(x)[
                                     contrast_name,
                                   ],
                                   p_value = get_p_values(x)[contrast_name, ],
                                   p = 1, q = 1, r = 0,
                                   cex = c(0.4, 1.5),
                                   col = c("#33333333", "#FF0000", "#FF666633"),
                                   pch = 19, feature_label = "feature",
                                   ylim = NULL,
                                   ...) {
  object <- x
  if (!(contrast_name %in% rownames(get_p_values(object)))) {
    stop(paste(
      "Choose a contrast in",
      paste(rownames(get_p_values(object)),
        collapse = ", "
      )
    ))
  }
  fold_change <- force(fold_change)
  p_value <- force(p_value)

  if (object$input$type == "1 sample") {
    stop("Can't do a volcano plot for one-sample tests!")
  }
  m <- object$input$n_dimensions
  stopifnot(m == length(fold_change))
  stopifnot(m == length(p_value))
  p_values_bound <- get_p_values(object)[contrast_name, ]
  thr <- thresholds(object)[1:m] # we select at most m hypotheses here

  vp <- volcano_plot(
    x = fold_change, p_values = p_value,
    thr = thr, p_values_bound = p_values_bound,
    p = p, q = q, r = r,
    cex = cex,
    col = col,
    pch = pch,
    feature_label = feature_label,
    ylim = ylim,
    add_signed_selections = TRUE#, ...
  )

  return(vp)
}
