#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-24 16:05:48
# @DESCRIPTION: filename
# @VERSION: v0.0.1

# Library -----------------------------------------------------------------

suppressPackageStartupMessages(library(magrittr))
library(ggplot2)
library(patchwork)
library(prismatic)
library(paletteer)
library(data.table)
#library(rlang)
library(glue)
library(parallel)
library(GetoptLong)
library(logger)
library(scales)
library(fs)

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

# header ------------------------------------------------------------------

# future: :plan(future: :multisession, workers = 10)

# load data ---------------------------------------------------------------

# load conn ---------------------------------------------------------------
conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)
# DBI::dbDisconnect(conn_all_hetero_af)
gseid_srrid_srrdir <- dplyr::tbl(
  conn_all_hetero_af,
  "gseid_srrid_srrdir"
) |>
  as.data.table()

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
fn_load_sc_and_sct <- function(.filepath, thegseid, thesrrid, ...) {
  .sc <- import(.filepath)
  sc_azimuth <- .sc$sc_azimuth
  rm(.sc)
  gc()
  sc_azimuth@meta.data |>
    tibble::rownames_to_column("barcode") |>
    as.data.table() |>
    # dplyr::left_join(
    #   forplot_ |>
    #     as.data.table() |>
    #     dplyr::mutate(
    #       barcode = as.character(barcode)
    #     ),
    # ) |>
    dplyr::mutate(
      barcode_new = glue::glue("{thegseid}-{thesrrid}-{barcode}")
    ) |>
    as.data.frame() -> d_merge

  new_names <- setNames(d_merge$barcode_new, d_merge$barcode)

  sc_azimuth <- RenameCells(
    sc_azimuth,
    new.names = new_names
  )

  sc_azimuth@meta.data <- d_merge |>
    tibble::column_to_rownames("barcode_new")

  sc_azimuth <- Seurat::SCTransform(
    sc_azimuth,
    assay = "RNA",
  )
  DefaultAssay(sc_azimuth) <- "SCT"

  sc_azimuth[["SCT"]]@scale.data <- matrix()
  sc_azimuth
}


fn_sct <- function(
  thegseid,
  thesrrid,
  ...
) {
  # thegseid <- gseid_srrid_srrdir$gseid[[1]]
  # thesrrid <- gseid_srrid_srrdir$srrid[[1]]
  library(Seurat)
  .dir <- path(
    "/home/liuc9/github/scMOCHA-data/data/",
    thegseid,
    "final",
    thesrrid
  )
  .dir_de <- path(
    .dir,
    "de"
  )

  dir_create(.dir_de)

  .sct_filepath <- path(
    .dir_de,
    "sc_azimuth.sct.qs"
  )

  sc_azimuth <- if (file_exists(.sct_filepath)) {
    log_fatal("{.sct_filepath} exists, skip!" |> glue::glue())
    file_delete(.sct_filepath)
    return(NULL)
  } else {
    sc_azimuth <- fn_load_sc_and_sct(
      .filepath = path(
        .dir,
        "sc_azimuth.rds.gz"
      ),
      thegseid = thegseid,
      thesrrid = thesrrid
    )
    export(
      sc_azimuth,
      file = .sct_filepath
    )
    sc_azimuth
  }
}


# body --------------------------------------------------------------------

gseid_srrid_srrdir |>
  dplyr::mutate(
    sct = parallel::mcmapply(
      FUN = fn_sct,
      thegseid = gseid,
      thesrrid = srrid,
      mc.cores = 20
    )
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
