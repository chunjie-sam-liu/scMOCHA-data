#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-06 10:58:26
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


# future: :plan(future: :multisession, workers = 10)


# load data ---------------------------------------------------------------
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/celltype-specific-variant"

gseid_srrid_ks_load <- import(
  file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test",
    "a_gseid_srrid_ks_load.qs"
  )
)

# function ----------------------------------------------------------------


# body --------------------------------------------------------------------

source("./analysis/00-colors.R")

gseid_srrid_ks_load |>
  dplyr::filter(p.value < 0.05) |>
  ggplot(
    aes(x = statistic)
  ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 100,
    fill = "grey50",
    color = "black",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    color = "red"
  )

gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 100
  )


thevariant <- "7833T>C"
thevariant <- "4175G>A"
thevariant <- "2101C>A"
thevariant <- "3664G>A"
thevariant <- "3243A>G"

gseid_srrid_ks_load |>
  dplyr::filter(variant == thevariant) |>
  dplyr::filter(
    p.value < 0.05,
    # statistic > 100
  ) |>
  dplyr::slice(1) |>
  tidyr::unnest(cols = celltype_af) |>
  dplyr::filter(af > 0) |>
  dplyr::mutate(
    celltype = gsub(
      "_",
      " ",
      celltype
    )
  ) |>
  dplyr::mutate(
    celltype = factor(celltype, levels = names(color_celltype) |> rev())
  ) |>
  ggplot(aes(
    x = af,
    y = celltype,
    fill = celltype
  )) +
  ggjoy::geom_joy(
    # scale = 0.9,
    # alpha = 0.8,
    rel_min_height = 0.01,
    size = 0.1
  ) +
  scale_fill_manual(
    values = color_celltype,
    na.value = "grey50"
  ) +
  ggjoy::theme_joy()


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
