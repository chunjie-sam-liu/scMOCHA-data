load_pkg(jutils)

# lintr::lint()
1 -> cj


lint(
  text = "foo <- function() { x <- 1
  1 -> cj }",
  linters = linters_with_defaults(
    commas_linter = commas_linter(),
    cyclocomp_linter = cyclocomp_linter(complexity_limit = 25),
    T_and_F_symbol_linter = T_and_F_symbol_linter(),
    duplicate_argument_linter = duplicate_argument_linter(),
    line_length_linter = NULL,
    infix_spaces_linter = NULL,
    spaces_left_parentheses_linter = NULL,
    spaces_inside_linter = NULL,
    indentation_linter = NULL,
    object_usage_linter = NULL,
    object_name_linter = NULL,
    assignment_linter = NULL,
    commented_code_linter = NULL,
    return_linter = NULL,
    brace_linter = NULL
  )
)
