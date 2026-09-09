#!/usr/bin/env Rscript
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-05-27
# @DESCRIPTION: Audit a named color vector from an R color file against expected levels.
# @VERSION: v0.1.0

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) {
    return(default)
  }
  args[[hit + 1]]
}

has_flag <- function(flag) {
  flag %in% args
}

if (has_flag("--help") || has_flag("-h")) {
  cat(
    "Usage: audit_color_file.R --file src/color.R --object celltype_colors [--levels A,B,C] [--output audit.pdf]\n"
  )
  quit(status = 0)
}

file <- get_arg("--file")
object <- get_arg("--object")
levels_arg <- get_arg("--levels", "")
output <- get_arg("--output")

if (is.null(file) || !file.exists(file)) {
  stop("Provide an existing --file.", call. = FALSE)
}

if (is.null(object) || !nzchar(object)) {
  stop("Provide --object with the named color vector to audit.", call. = FALSE)
}

env <- new.env(parent = globalenv())
sys.source(file, envir = env)

if (!exists(object, envir = env, inherits = FALSE)) {
  stop("Object not found in color file: ", object, call. = FALSE)
}

colors <- get(object, envir = env, inherits = FALSE)

if (is.function(colors)) {
  colors <- colors(10)
}

if (is.null(names(colors)) || any(!nzchar(names(colors)))) {
  cat("named_vector: false\n")
} else {
  cat("named_vector: true\n")
}

cat("object:", object, "\n")
cat("colors:", length(colors), "\n")

if (nzchar(levels_arg)) {
  expected <- strsplit(levels_arg, ",", fixed = TRUE)[[1]]
  expected <- trimws(expected)
  missing <- setdiff(expected, names(colors))
  extra <- setdiff(names(colors), expected)

  cat(
    "missing_levels:",
    if (length(missing)) paste(missing, collapse = ", ") else "none",
    "\n"
  )
  cat(
    "extra_colors:",
    if (length(extra)) paste(extra, collapse = ", ") else "none",
    "\n"
  )
}

if (!is.null(output) && nzchar(output)) {
  grDevices::pdf(output, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  if (requireNamespace("prismatic", quietly = TRUE)) {
    plot(prismatic::color(colors), main = object)
    plot(prismatic::clr_deutan(colors), main = paste(object, "deutan"))
    plot(prismatic::clr_protan(colors), main = paste(object, "protan"))
    plot(prismatic::clr_tritan(colors), main = paste(object, "tritan"))
  } else {
    plot.new()
    text(0.5, 0.5, "Install prismatic for color previews.")
  }
}
