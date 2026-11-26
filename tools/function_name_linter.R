function_name_linter <- function(styles = "snake_case") {

  #Define allowed styles
  regexes <- c(
    snake_case = "^[[:lower:]][[:lower:][:digit:]]*(_[[:lower:][:digit:]]+)*$",
    camelCase  = "^[[:lower:]][[:lower:][:digit:]]*([[:upper:]][[:lower:][:digit:]]*)*$",
    CamelCase  = "^[[:upper:]][[:lower:][:digit:]]*([[:upper:]][[:lower:][:digit:]]*)*$",
    SNAKE_CASE = "^[[:upper:]][[:upper:][:digit:]]*(_[[:upper:][:digit:]]+)*$",
    dotted.case = "^[[:lower:]][[:lower:][:digit:]]*(\\.[[:lower:][:digit:]]+)*$",
    lowercase  = "^[[:lower:]][[:lower:][:digit:]]*$",
    symbols    = "^[^[:alnum:]]+$"
  )

  #check style argument
  styles <- match.arg(styles, names(regexes), several.ok = TRUE)

  #use Linter function to create a Linter using
  lintr::Linter(
    function(source_expression) { #source_expression = code to evaluate
      #check lintr is used to file level
      if (!lintr::is_lint_level(source_expression, "file")) {
        return(list())
      }

      # find function names
      # nom <- function(...), nom = function(...), nom <<- \(x) ...
      xpath <- paste0(
        "//expr[LEFT_ASSIGN or EQ_ASSIGN]", #find affectation with <- or =
        "[expr[2][FUNCTION or OP-LAMBDA]]/expr[1]/SYMBOL" #filter affectation
        #where the second part is a function
      )
      #extract function names
      nodes <- xml2::xml_find_all(source_expression$full_xml_parsed_content, xpath)

      # on retire backticks et points initiaux (fonctions "privées" : .ma_fonction)
      noms <- gsub("^`|`$", "", xml2::xml_text(nodes))
      noms <- sub("^\\.+", "", noms)

      #check if function name corresponds to styles
      conforme <- vapply(
        noms,
        function(nom) any(vapply(regexes[styles], grepl, logical(1), x = nom)),
        logical(1)
      )

      #create lintr message
      lintr::xml_nodes_to_lints(
        nodes[!conforme],
        source_expression = source_expression,
        lint_message = paste0(
          "Function names must follow the style guidelines: ",
          paste(styles, collapse = ", "), "."
        ),
        type = "style"
      )
    },
    name = "function_name_linter",
    linter_level = "file"
  )
}

# function_name_linter = (function() {
#   source("tools/function_name_linter.R", local = TRUE)
#   function_name_linter(styles = c("snake_case", "symbols"))
# })(),
