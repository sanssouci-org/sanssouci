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
