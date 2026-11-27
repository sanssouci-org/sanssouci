set_struct.sanssouci_data <- function(
  object,
  forest_structure = NULL,
  leaf_list = NULL,
  s = NULL,
  H = NULL,
  flavor = 2
) {
  m <- get_n_hyp(object)
  if (!is.null(forest_structure) && !is.null(leaf_list)) {
    object$struct <- list(
      leaf_list = leaf_list,
      forest_structure = forest_structure)
  } else if (!is.null(leaf_list)) {
    object$struct <- dyadic_from_leaf_list(leaf_list, flavor)
  } else if (!is.null(s)) {
    object$struct <- dyadic_from_window_size(m, s, flavor)
  } else {
    object$struct <- dyadic_from_height(m, H, flavor)
  }
  object
}

new_ref_family_structured <- function(
  struct,
  zetas_list,
  alpha,
  method,
  step_down = FALSE,
  ...
) {
  parameters <- list(
    alpha = alpha,
    setting = "structured",
    method = method,
    step_down = step_down
  )
  parameters <- c(parameters, list(...))
  obj <- list(
    parameters = parameters,
    struct = struct,
    zetas_list = zetas_list
    )
  obj <- structure(
    obj,
    class = "ref_family_structured"
  )
  obj
}

switch_method_name <- function(method) {
  switch(method,
         "DKWM" = structured_zeta_dkwm,
         "HB" = structured_zeta_hb,
         "kBonf" = structured_zeta_kbonf,
         "trivial" = structured_zeta_trivial
  )
}

method_to_function <- function(method){
  if (typeof(method) == "list") {
    FUN <- list()
    for (h in seq_along(method)) {
      if (typeof(method[[h]]) == "list") {
        FUN[[h]] <- lapply(method[[h]], switch_method_name)
      } else {
        FUN[[h]] <- switch_method_name(method[[h]])
      }
    }  
  } else {
    FUN <- switch_method_name(method)
  }
  FUN
}

fit_structured.sanssouci_data <- function(
  object,
  alpha,
  method,
  step_down = FALSE,
  ...
) {
  struct <- get_struct(object)
  leaf_list <- struct$leaf_list
  forest_structure <- struct$forest_structure
  p_values <- get_p_values(object)
  FUN <- method_to_function(method)
  zetas_list <- structured_zetas(
    forest_structure,
    leaf_list,
    FUN,
    p_values,
    alpha,
    step_down = step_down,
    verbose = FALSE,
    ...)
  new_ref_family_structured(
    struct,
    zetas_list,
    alpha,
    method,
    step_down,
    ...
  )
}
