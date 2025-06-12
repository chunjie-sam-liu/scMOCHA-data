#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-11 11:48:13
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
library(glue)
library(parallel)
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
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant_cell.duckdb.1.2.1"
)

# function ----------------------------------------------------------------


# body --------------------------------------------------------------------
dplyr::tbl(conn, "all_variant_cell") |>
  dplyr::filter(
    variant_in_cell_cluster == "cell"
  ) |>
  dplyr::select(
    gseid, srrid, variant
  ) |>
  dplyr::distinct() |>
  as.data.table() ->
gseid_srrid_variant

all_variant_cell_table <- dplyr::tbl(conn, "all_variant_cell")

gseid_srrid_variant |>
  # head(100) |>
  dplyr::mutate(
    co = parallel::mcmapply(
      .x = srrid,
      .y = variant,
      FUN = \(.x, .y) {
        # .x <- "GSM4762179"
        # .y <- "11251A>G"

        log_trace(
          glue::glue(
            "Processing variant {.y} for srrid {.x}"
          )
        )
        all_variant_cell_table |>
          dplyr::filter(
            variant == .y,
            srrid == .x,
            variant_in_cell_cluster == "cell"
          ) |>
          dplyr::select(
            barcode, af, depth, variant_type, celltype
          ) |>
          as.data.table() ->
        .d
        .d |>
          dplyr::group_by(celltype) |>
          dplyr::summarise(sum_depth = sum(depth, na.rm = TRUE), mean_depth = mean(depth, na.rm = T)) ->
        .dd
        log_trace("has data in database ", nrow(.d))
        .d |>
          dplyr::count(
            celltype, variant_type
          ) |>
          dplyr::left_join(
            .dd,
            by = "celltype"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) ->
gseid_srrid_variant_co


gseid_srrid_variant_co |>
  tidyr::unnest(cols = co) |>
  tidyr::pivot_wider(
    names_from = variant_type,
    values_from = c(n),
  ) |>
  tidyr::nest(
    .by = c(gseid, srrid, variant),
    .key = "variant_celltype"
  ) ->
gseid_srrid_variant_celltype



gseid_srrid_variant_celltype |>
  dplyr::mutate(
    n_colorful = parallel::mcmapply(
      .x = variant_celltype,
      FUN = \(.x) {
        .x |>
          tidyr::pivot_longer(
            cols = -c(celltype, sum_depth, mean_depth),
            names_to = "group",
            values_to = "n",
          ) |>
          dplyr::mutate(
            n = ifelse(
              n >= 4,
              n,
              NA_real_
            )
          ) |>
          dplyr::filter(
            !is.na(n)
          ) |>
          dplyr::count(group) |>
          tidyr::pivot_wider(
            names_from = group,
            values_from = n,
            names_prefix = "n_"
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) |>
  tidyr::unnest(n_colorful) ->
gseid_srrid_variant_celltype_n


gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    !is.na(n_black),
    n_black == 8,
    n_colorful < 2
  ) |>
  # dplyr::slice(6) |>
  dplyr::filter(
    srrid == "GSM7080031"
  ) |>
  tidyr::unnest(cols = variant_celltype)



# ? real somatic mutation --------------------------------------------------------------------

ALLVARIANTS <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "all_variant.rds"
)) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  )




gseid_srrid_variant_celltype_n |>
  dplyr::filter(
    # srrid == "GSM7080031"
    variant %in% ALLVARIANTS$variant
  ) |>
  dplyr::filter(
    n_black >= 6,
    n_colorful < 6
  ) |>
  dplyr::filter(
    # srrid == "GSM7080027"
  ) |>
  dplyr::group_by(srrid) |>
  dplyr::filter(dplyr::n() > 3) |>
  dplyr::ungroup() |>
  dplyr::slice(3) |>
  tidyr::unnest(cols = variant_celltype)



# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
