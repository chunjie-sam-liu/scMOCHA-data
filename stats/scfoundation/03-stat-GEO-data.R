#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-26 15:37:33
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
project_source_sra <- readr::read_rds("/home/liuc9/github/scMOCHA-data/data/scfoundation/out/project_source_sra.rds.gz")

project_source_sra |> dplyr::glimpse()

# body --------------------------------------------------------------------
project_source_sra |>
  dplyr::filter(proj_source == "GEO") |>
  dplyr::select(proj_ID, source_name) |>
  dplyr::distinct() ->
project_source_sra_proj_ID_source_name


project_source_sra_proj_ID_source_name |>
  dplyr::count(source_name) |>
  dplyr::arrange(-n) |>
  dplyr::filter(
    grepl(
      pattern = "pbmc|Periperhal blood",
      x = source_name,
      ignore.case = TRUE
    )
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
