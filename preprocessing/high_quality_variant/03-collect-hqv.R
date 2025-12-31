#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-02 09:54:49
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
SRR <- data.table::fread(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/gse_srrid_srrdir.csv"
)
basedir <- "/home/liuc9/github/scMOCHA-data/data"
hqvdir <- file.path(basedir, "high_quality_variant")

# body --------------------------------------------------------------------
SRR |>
  dplyr::mutate(
    srrdir_hqv = file.path(
      hqvdir,
      gseid,
      "final",
      srrid
    )
  ) -> SRR_hqv

data.table::fwrite(
  SRR_hqv,
  file.path(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data",
    "gse_srrid_srrdir_hqv.csv"
  )
)

SRR_hqv |>
  dplyr::mutate(
    load_hqv = purrr::map(
      srrdir_hqv,
      ~ {
        data.table::fread(
          file.path(.x, "high_quality_variant.tsv")
        )
      }
    )
  ) -> SRR_hqv_load

readr::write_rds(
  SRR_hqv_load,
  file.path(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data",
    "gse_srrid_srrdir_hqv_load.rds"
  )
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
