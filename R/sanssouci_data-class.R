#' Create an object of class 'sanssouci_data'
#'
#' @param Y A matrix of \eqn{D} variables (hypotheses) by \eqn{n} observations
#' @param X A design matrix of \eqn{n} observations by \eqn{p} explanatory
#'   variables in the linear model
#' @param contrast A contrast matrix of size \eqn{L} tested contrasts (in row)
#'   and \eqn{p} columns corresponding to the parameters to be tested
#' @param truth An optional numeric vector of \eqn{m} values in \eqn{0, 1},
#'   the status of each null hypothesis (0 means H0 is true, 1 means H1 is
#'   true). Typically used in simulations.
#' @return An object of class "SansSouci"
#' @export
#' @examples
#' data(expr_ALL, package = "sanssouci.data")
#' groups <- ifelse(colnames(expr_ALL) == "NEG", 0, 1)
#' table(groups)
#' a <- new_sanssouci_data(Y = expr_ALL, X = groups)
#'
#' b <- new_sanssouci_data(Y = expr_ALL, X = groups, test = row_wilcoxon_tests)
#'
new_sanssouci_data <- function(Y = NULL,
                               p_values = NULL,
                               X = NULL,
                               contrast = NULL,
                               struct = NULL,
                               test_fun = NULL,
                               alternative = c("two.sided", "less", "greater"),
                               truth = NULL) {
  # first check: not all NULL
  if (is.null(Y) & is.null(p_values) & is.null(X) & is.null(contrast) &
      is.null(struct)) {
    stop(
      "At least one of 'Y', 'p_values', 'X', 'contrast', 'struct' ",
      "should be provided"
    )
  }
  alternative <- match.arg(alternative)
  if (is.null(struct)) {
    # Are the arguments present?
    signature <- c(!is.null(X), !is.null(contrast))

    n_obs <- ncol(Y)
    D <- nrow(Y)

    if (is.null(Y)) {
      ## can only do Simes and if so, need p-values
      ## setting = kfwer, method = Simes
      ## setting kfwer, method = lambdacalibration and n_permutations = 0
      if (is.null(p_values)) {
        stop("'p_values' should be provided if 'Y' is NULL")
      }
      if (!is.null(X) || !is.null(contrast) || !is.null(test_fun)) {
        message("If 'Y' is not provided, 'X', 'contrast' and 'test_fun'
                cannot be used")
      }
    } else { ## Y is not null

      if (is.null(X)) {
        X <- rep(1, n_obs)
        message("No design matrix provided: one-sample tests will be performed")
        # one sample test + message
        ## check consistency of function row tests ?
      }

      ## Y and X are not null

      if (is.null(contrast)) {
        ## No 'contrast': detection of one, two samples or continuous

        ## assume that if contrast is NULL, then X must be a vector.
        if(is.matrix(X)){
          stop("If contrast is NULL, please provide a vector for X.",
               "To test a matrix X, please provide a contrast matrix.")
        }

        n_groups <- length(unique(X))
        if (n_groups > 1) {
          .check_groups(X, ncol(Y))
        }


        fun_name <- as.character(substitute(test_fun))

        if (n_groups == 1) {
          type <- "1 sample"
          contrast_name <- "One sample test"
          if(is.null(test_fun)){
            test_fun <- row_z_tests
            fun_name <- "row_z_tests"
          }
        } else if (n_groups == 2) {
          type <- "2 samples"
          contrast_name <- "Group 1 vs group 0"
          if(is.null(test_fun)){
            test_fun <- row_welch_tests
            fun_name <- "row_welch_tests"
          }
          # } else if (n_groups == n) {
        } else { # by default consider a continuous
          type <- "Continuous case, correlation test"
          contrast_name <- "Correlation test"
          if (n_groups != n_obs) {
            message(
              "Ties were detected in the grouping variable, but a continuous",
              " vector was provided for groups, and a correlation test is being",
              " performed. If you intend to test associations across more than",
              " two groups, please provide a full design matrix using the X",
              " argument instead."
            )
          }
          groups <- X
          if(is.null(test_fun)){
            test_fun <- row_pearson_correlation_tests
            fun_name <- "row_pearson_correlation_tests"
          }
        }
        L <- 1
        P <- 1
      } else { # 'Y', 'X' and 'contrast' are not NULL: linear model

        .check_lm_test(t(Y), X, contrast)
        type <- "linear model"
        L <- nrow(contrast)
        P <- ncol(X)

        contrast_name <- rownames(contrast)
        if (is.null(contrast_name)) {
          contrast_name <- paste("contrast", seq_len(L))
          rownames(contrast) <- contrast_name
        }

        row_lm_test_with_x_c <- function(...) {
          row_lm_test(..., X = X, C = C)
        }
        test_fun <- row_lm_test_with_x_c
        fun_name <- "row_lm_test"
        groups <- matrix(1:n, ncol = 1)

      }

      if (!is.null(p_values)) {
        warning("p_values is provided as well as Y, then p_values will be",
                "discarded.")
      }

      res_test <- test_fun(Y, groups, alternative = alternative)
      p_values <- res_test$p_value
      estimates <- res_test$estimate #foldchange info


    }
  } else {
    # check struct format
    if (is.null(p_values)) {
      stop("Argument 'p_values' cannot be NULL if 'struct' is not NULL")
    }
    if (!is.null(Y) | !is.null(X) | !is.null(contrast) | !is.null(test_fun)) {
      warning(
        "Arguments 'Y', 'X', 'contrast', 'test_fun' not used ",
        "if 'struct' is not NULL"
      )
    }
  }
  if (!is.null(truth)) {
    all_0 <- identical(truth, rep(0, nrow(Y)))
    all_1 <- identical(truth, rep(1, nrow(Y)))
    if (!all_0 && !all_1) {
      .check_groups(truth, nrow(Y))
    }
  }
  obj <- structure(
    list(
      Y = Y,
      X = X,
      p_values = p_values,
      estimates = estimates,
      groups = groups,
      C = contrast,
      contrast_name = contrast_name,
      struct = struct,
      n_dimensions = D,
      n_contrasts = L,
      n_hypotheses = D * L,
      n_obs = n_obs,
      n_variables = P,
      test_fun = test_fun,
      test_fun_name = fun_name,
      type = type,
      truth = truth
    ),
    class = "sanssouci_data"
  )
  obj
}

get_p_values <- function(object) UseMethod("get_p_values")
#' `get_p_values`: get p-values
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_p_values.sanssouci_data <- function(object) {
  object$p_values
}

get_struct <- function(object) UseMethod("get_struct")
#' `get_struct`: get struct
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_struct.sanssouci_data <- function(object) {
  object$struct
}

get_n_hyp <- function(object) UseMethod("get_n_hyp")
#' `get_n_hyp`: get number of hypotheses
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_n_hyp.sanssouci_data <- function(object) {
  object$n_hyp
}

get_n_obs <- function(object) UseMethod("get_n_obs")
#' `get_n_obs`: get number of observations
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_n_obs.sanssouci_data <- function(object) {
  object$n_obs
}

get_test_fun <- function(object) UseMethod("get_test_fun")
#' `get_test_fun`: get test function
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_test_fun.sanssouci_data <- function(object) {
  object$test_fun
}

get_test_type <- function(object) UseMethod("get_test_type")
#' `get_test_type`: get type of test (one sample, two sample, correlation test,
#' linear model test)
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_test_type.sanssouci_data <- function(object) {
  object$type
}

get_groups_to_permute <- function(object) UseMethod("get_groups_to_permute")
#' `get_groups_to_permute`: get groups
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_groups_to_permute.sanssouci_data <- function(object) {
  object$groups
}

get_y <- function(object) UseMethod("get_y")
#' `get_y`: get Y
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_y.sanssouci_data <- function(object) {
  object$Y
}

get_truth <- function(object) UseMethod("get_truth")
#' `get_truth`: get truth
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_truth.sanssouci_data <- function(object) {
  object$truth
}



#' Print 'sanssouci_data' objects
#'
#' @rdname sanssouci_data-methods
#' @param x An object of class `sanssouci_data`
#' @param ... Not used
#' @param verbose Should detailed output be printed? Defaults to FALSE
#' @importFrom utils str
#' @export
print.sanssouci_data <- function(x, ..., verbose = FALSE) {
  object <- x
  rm(x)
  cat("'sanssouci_data' object:\n")
  cat("\tNumber of hypotheses:\t", get_n_hyp(object), "\n")
  if (!is.null(get_n_obs(object))) {
    cat("\tNumber of observations:\t", get_n_obs(object), "\n")
  }
  cat("\t", get_test_type(object), "\n", sep = "")
  cat("\tTest function:\t\t", object$test_fun_name, "\n", sep = "")
  if (verbose) {
    cat("Data:")
    cat("\n")
    str(object$Y)
    cat("\n")
    str(object$X)
    cat("\n")
    str(object$contrast)
    cat("\n")
  }
  truth <- object$truth
  if (!is.null(truth)) {
    cat("Truth: ")
    cat("\n")
    cat("\t", sum(truth), "false null hypotheses (signals)")
    pi0 <- 1 - mean(truth)
    cat(" out of ", get_n_hyp(object), " (pi0=", round(pi0, 3), ")", sep = "")
    cat("\n")
  }
invisible(object)
}

#' Fit Joint Error Rate Controlling Reference Family
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @param alpha A numeric value in `[0,1]`, the target (JER) risk
#' @param setting A character value, the considered setting for JER control. 
#'   * `"kFWER"`: joint k-Family-Wise Error Rate control (Blanchard et al. 2020)
#'   * `"level sets"`: not implemented yet.
#'   * `"structured"`: structured hypotheses (Durand et al. 2020)
#' @param method A character value, the method used to obtain JER control for a 
#' given setting. See [fit_kfwer.sanssouci_data], [fit_levelsets.sanssouci_data] 
#' and [fit_structured.sanssouci_data].
#' @param adapt A boolean value. Not implemented yet.
#' @param step_down A boolean value: should the step-down version of the 
#'   reference family be computed? Defaults to `FALSE`
#' @references Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc confidence bounds on false positives using reference families.
#' Annals of Statistics, 48(3), 1281-1303.
#' @references Durand, G., Blanchard, G., Neuvial, P., & Roquain, E. (2020).
#' Post hoc false positive control for structured hypotheses.
#' @references Enjalbert-Courrech, N. & Neuvial, P. (2022). Powerful and
#' interpretable control of false discoveries in two-group differential
#' expression studies. Bioinformatics. doi: 10.1093/bioinformatcs/btac693
#' @export
fit.sanssouci_data <- function(object, alpha,
                               setting = c("kFWER", "level sets", "structured"),
                               method, adapt = FALSE,
                               step_down = FALSE, ...) {
  
  if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' should be a numeric value in [0,1]")
  }
  
  setting <- match.arg(setting)
  if (setting == "kFWER") {
    ref_fam <- fit_kfwer(object, alpha = alpha, method = method,
                         adapt = adapt, step_down = step_down, ...)
  } else if (setting == "structured") {
    ref_fam <- fit_structured(object, alpha = alpha, method = method,
                          adapt = adapt, step_down = step_down, ...)
    
  } else if (setting == "level sets"){
    ref_fam <- fit_structured(object, alpha = alpha, method = method,
                              adapt = adapt, step_down = step_down, ...)
  } 
  return(ref_fam)
  
}
