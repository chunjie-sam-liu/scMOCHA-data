#!/usr/bin/env Rscript
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-05-27
# @DESCRIPTION: Inspect installed paletteer palettes by type and optional search text.
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
    "Usage: inspect_paletteer.R [--type all|discrete|dynamic|continuous] [--search TEXT]\n"
  )
  quit(status = 0)
}

if (!requireNamespace("paletteer", quietly = TRUE)) {
  stop("The paletteer package is not installed.", call. = FALSE)
}

type <- get_arg("--type", "all")
search <- get_arg("--search", "")

collect_names <- function(kind) {
  data <- switch(
    kind,
    discrete = paletteer::palettes_d_names,
    dynamic = paletteer::palettes_dynamic_names,
    continuous = paletteer::palettes_c_names,
    stop("Unknown palette type: ", kind, call. = FALSE)
  )

  names <- if ("palette" %in% names(data) && "package" %in% names(data)) {
    paste(data$package, data$palette, sep = "::")
  } else {
    as.character(data[[1]])
  }

  if (nzchar(search)) {
    names <- names[grepl(search, names, ignore.case = TRUE)]
  }

  names
}

kinds <- if (identical(type, "all")) {
  c("discrete", "dynamic", "continuous")
} else {
  type
}

for (kind in kinds) {
  names <- collect_names(kind)
  cat("\n#", kind, "palettes (", length(names), ")\n", sep = "")
  cat(paste(names, collapse = "\n"))
  cat("\n")
}
