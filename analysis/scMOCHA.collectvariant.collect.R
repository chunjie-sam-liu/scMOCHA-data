#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-10 17:24:07
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
srr |>
  dplyr::mutate(
    load = parallel::mclapply(
      srrdir,
      function(srrdir) {
        import(
          file.path(
            srrdir,
            "variant_info_from_heatmap.qs"
          )
        )
      },
      mc.cores = 20
    )
  ) ->
srr_load

srr_load |>
  dplyr::select(-srrdir) |>
  tidyr::unnest(cols = load) ->
srr_load_unnest

export(
  srr_load_unnest,
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.csv",
  format = "both",
)
export(
  srr_load_unnest |> as.data.table(),
  file = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.qs",
)


srr_load_unnest <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.qs")

v <- packageVersion("duckdb")
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.{v}" |> glue::glue()
)
DBI::dbWriteTable(
  conn,
  "all_variant_cell",
  srr_load_unnest,
  temporary = FALSE,
  overwrite = TRUE
)
DBI::dbDisconnect(conn, shutdown = TRUE)




# ? don't run below --------------------------------------------------------------------

conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)
dplyr::tbl(conn, "all_variant_cell") |>
  dplyr::filter(
    variant == "73A>G",
    srrid == "GSM4712885",
    variant_in_cell_cluster == "cell"
  ) |>
  as.data.table() |>
  dplyr::count(
    celltype, variant_type
  ) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = n
  )
DBI::dbDisconnect(conn, shutdown = TRUE)
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)
