new_ref_family_structured <- function(
  struct,
  zetas_list,
  alpha,
  method,
  stepdown = FALSE,
  ...
) {
  parameters <- list(
    alpha = alpha,
    setting = "structured",
    method = method,
    stepdown = stepdown,
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
         "dkwm" = structured_zeta_dkwm,
         "hb" = structured_zeta_hb,
         "kbonf" = structured_zeta_kbonf,
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
  stepdown = FALSE,
  ...
) {
  struct <- get_struct(object)
  leaf_list <- struct$leaf_list
  forest_structure <- struct$forest_structure
  p_values <- get_pvalues(object)
  FUN <- method_to_function(method)
  zetas_list <- structured_zetas(
    forest_structure,
    leaf_list,
    FUN,
    p_values,
    alpha,
    stepdown = stepdown,
    verbose = FALSE,
    ...)
  new_ref_family_structured(
    struct,
    zetas_list,
    alpha,
    method,
    stepdown,
    ...
  )
}
