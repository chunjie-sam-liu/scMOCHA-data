#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Thu Aug 22 12:34:44 2024
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
gseid <- "GSE226602"
datadir <- file.path(
  basedir,
  gseid
)

srrid_list <- readr::read_lines(
  file.path(
    datadir,
    "{gseid}.srrid.list" |> glue::glue()
  )
)

pheno <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.pheno.csv" |> glue::glue()
  )
) |>
  dplyr::filter(
    geo_accession %in% srrid_list
  )


pheno |>
  dplyr::glimpse()

pheno |>
  dplyr::select(
    srrid = geo_accession,
    age = `age:ch1`,
    disease = `disease state:ch1`,
    genotype = `genotype:ch1`,
    gender = `Sex:ch1`
  ) |>
  dplyr::arrange(disease, age, gender) -> pheno_sel

data.table::fwrite(
  x = pheno_sel,
  file = file.path(
    datadir,
    "{gseid}.pheno.select.csv" |> glue::glue()
  )
)

# body --------------------------------------------------------------------

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
