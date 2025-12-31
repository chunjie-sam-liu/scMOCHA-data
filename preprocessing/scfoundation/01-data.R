#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-21 13:42:29
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
# library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]
Options:

<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------

# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------

fn_plot_pie <- function(.d, .colors = NULL) {
  .d |>
    dplyr::select(group = 1, n) |>
    dplyr::arrange(-n) |>
    dplyr::mutate(csum = rev(cumsum(rev(n)))) %>%
    dplyr::mutate(pos = n / 2 + dplyr::lead(csum, 1)) %>%
    dplyr::mutate(pos = dplyr::if_else(is.na(pos), n / 2, pos)) %>%
    dplyr::mutate(percentage = n / sum(n)) |>
    dplyr::mutate(group = factor(group, levels = group)) -> .dd

  .scalefill <- if (is.null(.colors)) {
    ggsci::scale_fill_aaas(
      name = NULL
    )
  } else {
    scale_fill_manual(
      name = NULL,
      values = .colors
    )
  }
  .scalecolor <- if (is.null(.colors)) {
    ggsci::scale_color_aaas(
      name = NULL
    )
  } else {
    scale_color_manual(
      name = NULL,
      values = .colors
    )
  }

  .dd |>
    ggplot(aes(
      x = "",
      y = n,
    )) +
    geom_bar(
      aes(fill = group),
      stat = "identity",
      width = 1,
      color = "white",
      show.legend = FALSE
    ) +
    .scalefill +
    ggrepel::geom_label_repel(
      aes(
        y = pos,
        label = glue::glue("{group}\n{n} ({scales::percent(percentage)})"),
        color = group,
      ),
      size = 6,
      nudge_x = 1,
      nudge_y = 0,
      show.legend = FALSE,
      max.overlaps = Inf,
    ) +
    .scalecolor +
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
      ),
      # legend.position = "none"
    )
}

# load data ---------------------------------------------------------------
project_filename <- "/home/liuc9/github/scMOCHA-data/data/scfoundation/Projects.xlsx"

project <- readxl::read_xlsx(project_filename) |> as.data.table()

# body --------------------------------------------------------------------

project |>
  dplyr::mutate(
    project_source = parallel::mclapply(
      X = project_ID,
      FUN = function(.project_ID) {
        .source_ <- stringr::str_split(.project_ID, "-")[[1]][1]
        .ID <- paste0(
          stringr::str_split(.project_ID, "-")[[1]][-1],
          collapse = "-"
        )
        tibble::tibble(proj_source = .source_, proj_ID = .ID)
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = project_source) |>
  dplyr::mutate(
    sample_source = parallel::mclapply(
      X = sample_ID,
      FUN = function(.sample_ID) {
        .source_ <- stringr::str_split(.sample_ID, "-")[[1]][1]
        .ID <- paste0(
          stringr::str_split(.sample_ID, "-")[[1]][-1],
          collapse = "-"
        )
        tibble::tibble(samp_source = .source_, samp_ID = .ID)
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(cols = sample_source) -> project_source

project_source |>
  dplyr::count(proj_source) |>
  fn_plot_pie() -> project_source_pie


ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_pie.pdf",
  plot = project_source_pie,
  width = 10,
  height = 10
)

project_source |>
  data.table::fwrite(
    file = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source.csv",
    sep = ",",
    quote = FALSE,
    row.names = FALSE
  )


project_source |>
  dplyr::select(proj_source, proj_ID) |>
  dplyr::distinct() |>
  dplyr::count(proj_source) |>
  fn_plot_pie() -> nrow(project_source_proj_ID_pie)

ggsave(
  filename = "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_proj_ID_pie.pdf",
  plot = project_source_proj_ID_pie,
  width = 10,
  height = 10
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
