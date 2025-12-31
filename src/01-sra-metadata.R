#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Aug 16 17:45:09 2024
# @DESCRIPTION: filename

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

# gseid <- "GSE163668"

basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
# basedir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"
verbose <- FALSE
spec <- "
Usage: Rscript foorbar.R [options]

Options:
<gseid=s> gseid
<basedir=s> basedir
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

fn_get_html <- function(theid, prefix, datadir) {
  # .url <- glue::glue("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={theid}") # geo
  # .url <- glue::glue("https://www.ncbi.nlm.nih.gov/Traces/study/?acc={theid}") # sraruntable
  .url <- glue::glue("{prefix}={theid}")

  # .html <- rvest::read_html(.url)
  httr::GET(.url)

  .of <- file.path(
    datadir,
    "{theid}.html" |> glue::glue()
  )

  xml2::write_html(
    x = .html,
    file = .of
  )
  .of
}

fn_parse_bioprojectid <- function(.of) {
  rvest::read_html(.of) |>
    rvest::html_elements("td") %>%
    .[rvest::html_text2(.) == "BioProject"] -> .ele

  .ele |>
    rvest::html_element(xpath = "./parent::tr") |>
    rvest::html_elements("td") %>%
    .[[2]] -> .the_ele

  .projid <- rvest::html_text2(.the_ele)
  .projid
}

fn_parse_sratable <- function(.of) {
  rvest::read_html(.of) |>
    rvest::html_table()
}

fn_edirect_gds_gseid <- function(.gseid, datadir = datadir) {
  .cmd <- "esearch -db gds -query {.gseid}|efetch -format docsum >{datadir}/{.gseid}.edirect.gds.xml" |>
    glue::glue()
  # system(.cmd)
  .cmd
}

# load data ---------------------------------------------------------------

# body --------------------------------------------------------------------

log_warn(gseid)
datadir <- file.path(
  basedir,
  gseid
)

dir.create(
  path = datadir,
  showWarnings = F,
  recursive = T
)

# get gsm -----------------------------------------------------------------

gse <- GEOquery::getGEO(
  GEO = gseid,
  destdir = datadir,
  GSEMatrix = TRUE
)


gse[[1]] |>
  Biobase::phenoData() |>
  Biobase::pData() -> gse_pheno

data.table::fwrite(
  x = gse_pheno,
  file = file.path(
    datadir,
    "{gseid}.pheno.csv" |> glue::glue()
  )
)

# cmd_edirect_gds_gseid <- fn_edirect_gds_gseid(gseid, datadir = datadir)
# readr::write_lines(
#   x = cmd_edirect_gds_gseid,
#   file = file.path(
#     datadir,
#     "00.edirect.gds.{gseid}.sh" |> glue::glue()
#   )
# )

# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
