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


# load data ---------------------------------------------------------------
project_filename <- "/home/liuc9/github/scMOCHA-data/data/scfoundation/Projects.xlsx"

project <- readxl::read_xlsx(project_filename) |> as.data.table()

# body --------------------------------------------------------------------


project |>
  tidyr::separate(col = project_ID, into = c("source", "ID"), sep = "-") |>
  dplyr::filter(source == "GEO") |>
  as.data.table() |>
  dplyr::select(ID) |>
  dplyr::distinct() ->
project_geo



# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
