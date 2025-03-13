#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-13 11:29:32
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


basedir <- "/home/liuc9/github/scMOCHA-data/data"


gseids = (
  gseids <- c(
    "GSE143353", # done
    "GSE147794", # finished, some are not success, need rerun
    "GSE148215", # under run, done
    "GSE153421", # under run, done
    "GSE163314", # under run, done
    "GSE163633", # under run, some are not success, need rerun
    "GSE164690", # under run, some are not success, need rerun
    "GSE165087", # not run
    "GSE165496", # under run, not run.
    "GSE165822", # under run, not run.
    "GSE167825", # under run, some are not success, need rerun
    "GSE168453", # under run, some are not success, need rerun
    "GSE174125", # under run, done
    "GSE184703" # under run, some are not success, need rerun
  )
)


gseids_meta_scfoundation <- tibble::tibble(
  GSE_ID = c(
    # scfoundation2
    "GSE143353", # done
    "GSE147794", # finished, some are not success, need rerun
    "GSE148215", # under run, done
    "GSE153421", # under run, done
    "GSE163314", # under run, done
    "GSE163633", # under run, some are not success, need rerun
    "GSE164690", # under run, some are not success, need rerun
    "GSE165087", # not run
    "GSE165496", # under run, not run.
    "GSE165822", # under run, not run.
    "GSE167825", # under run, some are not success, need rerun
    "GSE168453", # under run, some are not success, need rerun
    "GSE174125", # under run, done
    "GSE184703" # under run, some are not success, need rerun
  ),
) |>
  dplyr::mutate(
    samples = purrr::map(
      GSE_ID,
      .f = \(.x) {
        basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"
        data.table::fread(
          file.path(basedir, .x, "out", glue::glue("{.x}.cell_ratio_and_variant_clean.csv"))
        ) ->
        .d
        tibble::tibble(
          samples = nrow(.d),
          Disease = "-",
          Source = "PBMC",
          Chemistry = unique(.d$Chemistry)[[1]],
          Publication = "-"
        )
      }
    ),
  ) |>
  tidyr::unnest(cols = samples)
# body --------------------------------------------------------------------


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
