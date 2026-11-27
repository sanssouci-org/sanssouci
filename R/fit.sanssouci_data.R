fit.sanssouc_data <- function(object, alpha,
                              setting = c("kfwer", "levelets", "struct"),
                              method, adapt = FALSE,
                              stepdown = FALSE, ...) {

  if (mode(alpha) != "numeric" || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' should be a numeric value in [0,1]")
  }

  match.arg(setting)

  if(setting == "kfwer"){

    object <- fit_kfwer(object, alpha = alpha, method = method,
                        adapt = adapt, stepdown = stepdown, ...)
  } else if (setting == "struct") {
    object <- fit_struct(object, alpha = alpha, method = method,
                        adapt = adapt, stepdown = stepdown, ...)

  } else if (setting == "levelsets"){
    stop("TODO : the setting levelsets is not currently available")
  } else {
    stop("Please give good setting name")
  }

  return(object)


}
