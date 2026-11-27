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
#' a <- SansSouci(Y = expr_ALL, groups = groups)
#'
#' res <- fit(a, B = 100, alpha = 0.1)
#' plot(res, xmax = 500)
#' res_beta <- fit(res, B = 100, alpha = 0.1, family = "Beta", K = 10)
#' plot(res_beta, xmax = 500)
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
  if (is.null(Y) & is.null(p_values) & is.null(X) & is.null(contrast) 
      & is.null(struct)) {
    stop("At least one of 'Y', 'p_values', 'X', 'contrast', 'struct' ",
         "should be provided")
  }
  alternative <- match.arg(alternative)
  if (is.null(struct)) {
    # Are the arguments present?
    signature <- c(!is.null(X), !is.null(contrast))
    
    n_obs <- ncol(Y)
    D <- nrow(Y)
    
    if (is.null(Y)) {
      ## can only do Simes and if so, need p-values
      if (is.null(p_values)) {
        stop("'p_values' should be provided if 'Y' is NULL")
      }
      if (!is.null(X) | !is.null(contrast) | !is.null(test_fun)) {
        message("If 'Y' is not provided, 'X', 'contrast' and 'test_fun' 
                cannot be used")
      }
    } else { ## Y is not null
      if (is.null(X)) {
        # one sample test + message
        type <- "1 sample"
        contrast_name <- "One sample test"
        message("No design matrix provided: one-sample tests will be performed")
        ## check consistency of function row tests ?
      } else { ## Y and X are not null
        if (is.null(contrast)) {
          ## No 'contrast': detection of one, two samples or continuous
          n_groups <- length(unique(X))
          if (n_groups > 1) {
            .check_groups(X, ncol(Y))
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
          }
          L <- 1                              
          P <- 1
        } else { # 'Y', 'X' and 'contrast' are not NULL: linear model
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
        }
      }
    }
  } else {
    # check struct format    
    if (is.null(p_values)) {
      stop("Argument 'p_values' cannot be NULL if 'struct' is not NULL")
    }
    if (!is.null(Y) | !is.null(X) | !is.null(contrast) | !is.null(test_fun)) {
      warning("Arguments 'Y', 'X', 'contrast', 'test_fun' not used ", 
              "if 'struct' is not NULL")
    }
  }
  if (!is.null(truth)) {
    all_0 <- identical(truth, rep(0, nrow(Y)))
    all_1 <- identical(truth, rep(1, nrow(Y)))
    if (!all_0 && !all_1) {
      .check_groups(truth, nrow(Y))
    }
  }
  obj <- struct(
    list(
      Y = Y,
      X = X,
      p_values = p_values,
      groups = groups,
      C = contrast,
      contrast_name = contrast_name,
      struct = struct,
      n_dimensions = D,
      n_contrasts = L,
      n_hypotheses = D * L,
      n_obs = n,
      n_variables = P,
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
