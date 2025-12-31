#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-28 12:13:08
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

# load data ---------------------------------------------------------------

datadir <- "/home/liuc9/github/scMOCHA-data/data"
gseid <- "GSE181279"
gsm <- "GSM5494116"

# body --------------------------------------------------------------------

stats <- data.table::fread(
  file.path(datadir, gseid, "final", gsm, "cell.variant_stats.tsv.gz")
)


stats |>
  dplyr::filter(
    n_cells_conf_detected > 3
  )

variant_annotation <- data.table::fread(
  file.path(datadir, gseid, "final", gsm, "variant_annotation.tsv")
) |>
  dplyr::mutate(
    variant = glue::glue("{Position}{Ref}>{Alt}")
  )

stats |>
  dplyr::mutate(
    vmr_log = log10(vmr)
  ) |>
  dplyr::filter(
    strand_correlation > 0.3 & vmr_log > log10(0.01)
  ) -> stats_filtered


ggvenn::ggvenn(
  data = list(
    "variant_annotation" = variant_annotation$variant,
    "stats" = stats$variant,
    "stats_filtered" = stats_filtered$variant
  ),
)

#
#
#

stats |>
  dplyr::mutate(
    vmr_log = log10(vmr)
  ) |>
  dplyr::mutate(
    color = ifelse(
      strand_correlation > 0.65 &
        vmr_log > log10(0.01) &
        variant %in% variant_annotation$variant,
      "red",
      "black"
    )
  ) |>
  ggplot(aes(x = strand_correlation, y = vmr_log)) +
  geom_point(aes(color = color)) +
  scale_color_identity() +
  geom_vline(
    xintercept = 0.65,
    linetype = 20,
    color = "red"
  ) +
  geom_hline(
    yintercept = log10(0.01),
    linetype = 20,
    color = "red"
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
