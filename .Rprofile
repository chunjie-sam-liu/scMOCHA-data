# Blogdown options --------------------------------------------------------
options(
  blogdown.author = "Chun-Jie Liu",
  servr.daemon = FALSE,
  blogdown.ext = ".Rmd",
  blogdown.subdir = "post",
  blogdown.yaml.empty = TRUE,

  # General options ---------------------------------------------------------
  repos = c(CRAN = "https://cloud.r-project.org"),
  prompt = "Chun-Jie>",
  digits = 4,
  show.signif.stars = FALSE,
  stringsAsFactors = FALSE,

  # ggplot2 v3 options ------------------------------------------------------
  ggplot2.continuous.color = "viridis",
  ggplot2.continuous.fill = "viridis",
  ggrepel.max.overlaps = Inf,
  future.globals.maxSize = 100 * 1024^3,
  warn = 1,
  warnPartialMatchArgs = FALSE,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE,
  showWarnCalls = TRUE,
  showErrorCalls = TRUE,
  datatable.print.class = TRUE,
  datatable.print.colnames = "top",
  pillar.subtle = FALSE,
  pillar.neg = FALSE,
  # readr options ------------------------------------------------------
  readr.num_columns = 0L,
  readr.show_progress = FALSE,
  readr.show_col_types = FALSE,
  dplyr.summarise.inform = FALSE,
  gargle_oauth_email = TRUE,
  gargle_oauth_cache = TRUE,
  devtools.install.args = c("--no-multiarch", "--no-test-load"),
  styler.cache_root = "styler-perm",
  testthat.default_check_reporter = "progress",
  languageserver.formatting_style = function(.options) {
    style = styler::tidyverse_style(indent_by = .options$tabSize)
    style$token$force_assignment_op = NULL
    style
  },
  # rio options ------------------------------------------------------

  rio.import.class = "data.table"
)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
options(lintr.linter_file = ".lintr")

# Radian ----------------------------------------------------------------

# Do not copy the whole configuration, just specify what you need!
# see https://pygments.org/styles
# for a list of supported color schemes, default scheme is "native"
options(radian.color_scheme = "native")

# # either  `"emacs"` (default) or `"vi"`.
# options(radian.editing_mode = "emacs")
# # enable various emacs bindings in vi insert mode
# options(radian.emacs_bindings_in_vi_insert_mode = FALSE)
# # show vi mode state when radian.editing_mode is `vi`
# options(radian.show_vi_mode_prompt = TRUE)
# options(radian.vi_mode_prompt = "\033[0;34m[{}]\033[0m ")

# indent continuation lines
# turn this off if you want to copy code without the extra indentation;
# but it leads to less elegent layout
options(radian.indent_lines = TRUE)

# auto match brackets and quotes
options(radian.auto_match = TRUE)

# enable the [prompt_toolkit](https://python-prompt-toolkit.readthedocs.io/en/master/index.html) [`auto_suggest` feature](https://python-prompt-toolkit.readthedocs.io/en/master/pages/asking_for_input.html#auto-suggestion)
# this option is experimental and is known to break python prompt, use it with caution
options(radian.auto_suggest = FALSE)

# highlight matching bracket
options(radian.highlight_matching_bracket = FALSE)

# auto indentation for new line and curly braces
options(radian.auto_indentation = TRUE)
options(radian.tab_size = 2)

# pop up completion while typing
options(radian.complete_while_typing = TRUE)
# the minimum length of prefix to trigger auto completions
options(radian.completion_prefix_length = 2)
# timeout in seconds to cancel completion if it takes too long
# set it to 0 to disable it
options(radian.completion_timeout = 0.05)
# add spaces around equals in function argument completion
options(radian.completion_adding_spaces_around_equals = TRUE)

# automatically adjust R buffer size based on terminal width
options(radian.auto_width = TRUE)

# insert new line between prompts
options(radian.insert_new_line = TRUE)

# max number of history records
options(radian.history_size = 20000)
# where the global history is stored, environmental variables will be expanded
# note that "~" is expanded to %USERPROFILE% or %HOME% in Windows
# options(radian.global_history_file = "~/.radian_history")
# the filename that local history is stored, this file would be used instead of
# `radian.global_history_file` if it exists in the current working directory
# options(radian.local_history_file = ".radian_history")
# when using history search (ctrl-r/ctrl-s in emacs mode), do not show duplicate results
# options(radian.history_search_no_duplicates = FALSE)
# ignore case in history search
# options(radian.history_search_ignore_case = FALSE)
# do not save debug browser commands such as `Q` in history
# options(radian.history_ignore_browser_commands = TRUE)

# custom prompt for different modes
options(radian.prompt = "\033[0;34mradian>\033[0m ")
options(radian.shell_prompt = "\033[0;31m#!>\033[0m ")
options(radian.browse_prompt = "\033[0;33mBrowse[{}]>\033[0m ")

# stderr color format
options(radian.stderr_format = "\033[0;31m{}\033[0m")

# enable reticulate prompt and trigger `~`
options(radian.enable_reticulate_prompt = TRUE)


options(vsc.str.max.level = 2)

if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
  if (requireNamespace("httpgd", quietly = TRUE)) {
    options(vsc.plot = FALSE)
    options(device = function(...) {
      httpgd::hgd(silent = TRUE)
      .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
    })
  }
}

# if (interactive() && Sys.getenv("RSTUDIO") == "") {
#   source(file.path(
#     Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"),
#     ".vscode-R",
#     "init.R"
#   ))
# }

# Functions ----------------------------------------------------------------
tryCatch(
  {
    func_file <- if (file.exists("/workspace/func.R")) {
      "/workspace/func.R"
    } else if (file.exists("~/github/renv/func.R")) {
      "~/github/renv/func.R"
    } else {
      file.path(getwd(), "func.R")
    }
    source(func_file, local = TRUE)
  },
  error = function(e) {
    message("Failed to source func.R: ", e$message)
  }
)

## ensure zip works for openxlsx2 / rmarkdown / quarto
if (nzchar(Sys.which("zip"))) {
  options(zip = Sys.which("zip"))
}

zipcmd <- Sys.which("zip")
if (nzchar(zipcmd)) {
  Sys.setenv(R_ZIPCMD = zipcmd)
}

.First <- function() {
  if (interactive()) {
    suppressMessages({
      load_pkg(jutils)
    })
    dotenv()
    suppressMessages({
      conflicted::conflicts_prefer(dplyr::filter, fs::path)
    })
    # Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "FALSE")

    logger::log_threshold(logger::TRACE)
    logger::log_layout(logger::layout_glue_colors)

    message("           ❤️Hello, Chun-Jie.\n")
  }
}


.Last <- function() {
  if (interactive()) {
    message("\n\n           👋Bye, Chun-Jie.\n\n")
  }
}
