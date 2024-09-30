#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------

# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: Fri Sep 27 11:32:05 2024
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

# future::plan(future::multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
gseidlist <- c(
  "GSE163668",
  "GSE149689",
  "GSE155223",
  "GSE155673",
  "GSE157344",
  "GSE166992",
  "GSE171555"
)

basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
# body --------------------------------------------------------------------
length(gseidlist)
gseid <- gseidlist[[7]]
datadir <- file.path(basedir, gseid)
outdir <- file.path(
  datadir, "out"
)
dir.create(outdir, showWarnings = F,recursive = T)

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
srarun <- data.table::fread(
  file.path(
    datadir,
    "{gseid}.SraRunTable" |> glue::glue()
  )
)
pheno |> dplyr::glimpse()

srarun |> dplyr::glimpse()
gseid
srarun |>
  dplyr::select(
    srrid = `Sample Name`,
    age = dplyr::contains("age"),
    gender = dplyr::contains("gender"),
    sex = dplyr::contains("sex")
  )

pheno |>         #dplyr::glimpse()
  dplyr::select(
    srrid = geo_accession,
    age = dplyr::contains("age"),
    gender = dplyr::contains("gender"),
    # race = `characteristics_ch1.6`,
    sex = dplyr::contains("sex"),
    disease = `disease state:ch1`
  ) |>
  dplyr::mutate(
    disease = gsub("", "", x = disease)
  ) ->
  suc;suc

data.table::fwrite(
  x = suc,
  file = file.path(
    outdir,
    "{gseid}.meta.csv" |> glue::glue()
  )
)


# footer ------------------------------------------------------------------

# future::plan(future::sequential)

# save image --------------------------------------------------------------
