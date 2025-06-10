#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-10 12:32:56
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
srr_filename <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
srr <- import(srr_filename)

# function ----------------------------------------------------------------


# body --------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)
srr |>
  # head(10) |>
  dplyr::mutate(
    load = parallel::mclapply(
      X = srrdir,
      FUN = \(.srrdir) {
        log_info(
          sprintf(
            "Loading %s",
            .srrdir
          )
        )
        .sc <- import(
          file.path(.srrdir, "sc_azimuth.rds.gz")
        )

        .sc$sc_azimuth@meta.data |>
          tibble::rownames_to_column("barcode") |>
          dplyr::select(
            -c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, Percent.Largest.Gene),
            -dplyr::contains("score"),
          ) |>
          data.table::as.data.table() ->
        .d

        .colnames <- gsub("\\.", "_", gsub("predicted.", "", colnames(.d)))
        colnames(.d) <- .colnames
        export(
          .d,
          file.path(.srrdir, "sc_azimuth_celltype.csv")
        )

        log_success(
          sprintf(
            "Load %s, %d barcodes",
            .srrdir,
            nrow(.d)
          )
        )
        rm(.sc)
        gc()
        return(.d)
      },
      mc.cores = 50
    )
  ) ->
srr_load

srr_load |>
  dplyr::select(-srrdir) |>
  tidyr::unnest(cols = c(load)) ->
srr_load_unnest
export(
  srr_load_unnest,
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype_detail.csv",
  format = "both"
)

export(
  srr_load_unnest,
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype_detail.qs"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
