#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Sun Aug 18 21:18:32 2024
# @DESCRIPTION: filename

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(GetoptLong)
library(logger)

# args --------------------------------------------------------------------


gseid <- "GSE163668"

# s: string, i: integer, f: float, !: boolean
# @: array
# %: hash
# default: default value specified here.
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<gseid=s> gseid
<verbose!> Print messages
"

GetoptLong.options(help_style = "two-column")
GetoptLong(spec, template_control = list(opt_width = 21))

# src ---------------------------------------------------------------------


# header ------------------------------------------------------------------
log_threshold(TRACE)
log_layout(layout_glue_colors)

# log_info('Starting the script...')
# log_debug('This is the second log line')
# log_trace('Note that the 2nd line is being placed right after the 1st one.')
# log_success('Doing pretty well so far!')
# log_warn('But beware, as some errors might come :/')
# log_error('This is a problem')
# log_debug('Note that getting an error is usually bad')
# log_error('This is another problem')
# log_fatal('The last problem')

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------


basedir <- "/home/liuc9/github/scMOCHA-data/data"
datadir <- file.path(
  basedir, gseid
)

dir.create(
  path = datadir,
  showWarnings = F,
  recursive = T
)

runfile <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.runfile.csv" |> glue::glue()
  )
)

# body --------------------------------------------------------------------
gsm <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.metadata.gsm.csv" |> glue::glue()
  )
)

gsm |>
  dplyr::rename(
    srrid = run_accession
  ) |>
  dplyr::inner_join(
    runfile, by = "srrid"
  ) |>
  dplyr::group_by(experiment_name) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    rename = purrr::map2(
      .x = experiment_name,
      .y = data,
      .f = \(.x, .y) {
        gsmdir <- file.path(
          datadir,
          .x
        )

        dir.create(
          gsmdir,
          showWarnings = F,
          recursive = T
        )

        .x
        .y |>
          dplyr::mutate(
            sym
          )

      }
    )
  )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
