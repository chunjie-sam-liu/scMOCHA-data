#!/usr/bin/env Rscript
# @AUTHOR: Chunjie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-05-27
# @DESCRIPTION: Preview a paletteer palette with prismatic or base R plotting.
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
    "Usage: preview_palette.R --palette package::palette [--type discrete|dynamic|continuous] [--n 10] [--output palette_preview.pdf]\n"
  )
  quit(status = 0)
}

if (!requireNamespace("paletteer", quietly = TRUE)) {
  stop("The paletteer package is not installed.", call. = FALSE)
}

palette <- get_arg("--palette")
type <- get_arg("--type", "discrete")
n <- as.integer(get_arg("--n", "10"))
output <- get_arg("--output")

if (is.null(palette) || !nzchar(palette)) {
  stop("Provide --palette as package::palette.", call. = FALSE)
}

colors <- switch(
  type,
  discrete = as.character(paletteer::paletteer_d(palette)),
  dynamic = as.character(paletteer::paletteer_dynamic(palette, n)),
  continuous = as.character(paletteer::paletteer_c(palette, n = n)),
  stop("Unknown --type. Use discrete, dynamic, or continuous.", call. = FALSE)
)

if (identical(type, "discrete") && length(colors) > n) {
  colors <- colors[seq_len(n)]
}

plot_colors <- function(colors, title) {
  if (requireNamespace("prismatic", quietly = TRUE)) {
    plot(prismatic::color(colors), main = title)
  } else {
    old_par <- par(mar = c(1, 1, 3, 1))
    on.exit(par(old_par), add = TRUE)
    plot.new()
    plot.window(xlim = c(0, length(colors)), ylim = c(0, 1))
    rect(
      seq_along(colors) - 1,
      0,
      seq_along(colors),
      1,
      col = colors,
      border = NA
    )
    title(main = title)
  }
}

if (!is.null(output) && nzchar(output)) {
  grDevices::pdf(output, width = 8, height = 2)
  on.exit(grDevices::dev.off(), add = TRUE)
}

plot_colors(colors, paste(type, palette))
cat(paste(colors, collapse = "\n"), "\n", sep = "")
