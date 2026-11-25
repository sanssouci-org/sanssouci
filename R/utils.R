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


categ_check <- function(categ, n) {
  name <- as.character(substitute(categ))
  if (length(categ) != n) {
    stop(name, " should be of length ", n, ", not", length(categ))
  }
  categ <- as.factor(categ)
  cats <- levels(categ)
  if (!identical(cats, c("0", "1")) && length(cats) <= 2) {
    stop(
      "'", name,
      "' should consist only of '0' and '1' or distinct continuous values."
    )
  }
}
