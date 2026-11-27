#' Match the singular and plural forms of a word with a quantity
#'
#' @param word a character, the word to adapt
#' @param n a numeric, the quantity
#'
#' @returns a character
#'
pluralize <- function(word, n) {
  word <- sub("s$", "", word)

  if (n == 1)
    return(word)

  if (grepl("[sxz]$", word) || grepl("(ch|sh)$", word))
    return(paste0(word, "es"))

  if (grepl("[^aeiou]y$", word))
    return(sub("y$", "ies", word))

  paste0(word, "s")
}


.check_groups <- function(groups, n) {
  name <- as.character(substitute(groups))
  if (length(groups) != n) {
    stop(name, " should be of length ", n, ", not", length(groups))
  }
  groups <- as.factor(groups)
  cats <- levels(groups)
  if (!identical(cats, c("0", "1")) && length(cats) <= 2) {
    stop(
      "'", name,
      "' should consist only of '0' and '1' or distinct continuous values."
    )
  }
}
