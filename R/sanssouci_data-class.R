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
                               structure = NULL,
                               test_fun = NULL,
                               truth = NULL) {
  # first check: not all NULL
  if (is.null(structure)) {
    # Are the arguments present?
    signature <- c(!is.null(X), !is.null(contrast))
    
    n <- ncol(Y)
    D <- nrow(Y)
    
    # one sample case
    if (all(X == 1)) {
      if (!is.null(contrast)) {
        warning("Argument 'contrast' not used if 'structure' is not NULL")
      }
      groups <- rep(1, ncol(Y))
      type <- "1 sample"
      L <- 1
      P <- 1
      contrast_name <- "one sample test"
    } else if (!is.null(X) & is.null(contrast)) {
      ## Only groups is given: detection of one, two samples or continuous
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
    
  } else {
    # check structure format    
    if (is.null(p_values)) {
      stop("Argument 'p_values' cannot be NULL if 'structure' is not NULL")
    }
    if (!is.null(Y) | !is.null(X) | !is.null(contrast) | !is.null(test_fun)) {
      warning("Arguments 'Y', 'X', 'contrast', 'test_fun' not used ", 
              "if 'structure' is not NULL")
    }
  }
  if (!is.null(truth)) {
    all_0 <- identical(truth, rep(0, nrow(Y)))
    all_1 <- identical(truth, rep(1, nrow(Y)))
    if (!all_0 && !all_1) {
      categ_check(truth, nrow(Y))
    }
  }
  obj <- structure(
    list(
      Y = Y,
      X = X,
      p_values = p_values,
      groups = groups,
      C = contrast,
      contrast_name = contrast_name,
      structure = structure,
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

get_structure <- function(object) UseMethod("get_structure")
#' `get_structure`: get structure
#'
#' @rdname sanssouci_data-methods
#' @param object An object of class `sanssouci_data`
#' @export
get_structure.sanssouci_data <- function(object) {
  object$structure
}

