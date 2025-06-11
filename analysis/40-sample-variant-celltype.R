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
            "Processing variant {thevariant} for srrid {thesrrid}"
          )
        )
        all_variant_cell_table |>
          dplyr::filter(
            variant == .y,
            srrid == .x,
            variant_in_cell_cluster == "cell"
          ) |>
          dplyr::select(
            barcode, af, variant_type, celltype
          ) |>
          as.data.table() ->
        .d
        log_trace("has data in database ", nrow(.d))
        .d |>
          dplyr::count(
            celltype, variant_type
          ) |>
          tidyr::pivot_wider(
            names_from = variant_type,
            values_from = n
          ) |>
          dplyr::mutate(
            realvariant = ifelse(
              colorful >= 3 & !is.na(colorful),
              "yes",
              "no"
            )
          )
      },
      mc.cores = 20,
      SIMPLIFY = FALSE
    )
  ) ->
gseid_srrid_variant_co


gseid_srrid_variant_co |>
  dplyr::mutate(
    n = parallel::mcmapply(
      .x = co,
      FUN = \(.x) {
        # .x <- gseid_srrid_variant_co$co[[1]]
        .x |>
          dplyr::filter(realvariant == "yes") |>
          nrow()
      },
      SIMPLIFY = TRUE,
      mc.cores = 20
    )
  ) ->
gseid_srrid_variant_co_n

gseid_srrid_variant_co_n |>
  dplyr::select(n) |>
  ggplot(aes(x = n)) +
  geom_histogram(
    binwidth = 1,
    fill = "steelblue",
    color = "black"
  ) +
  labs(
    title = "Number of real variants per srrid",
    x = "Number of real variants",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


gseid_srrid_variant_co_n |>
  dplyr::filter(n == 1) |>
  dplyr::slice(2) |>
  tidyr::unnest(cols = co)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
