#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-05 16:28:44
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
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

all_heteroplasmic_af <- import(
  file.path(cleandatadir, "all_hetero_af.cluster.fst"),
) |>
  dplyr::select(-num_variants) |>
  dplyr::rename(celltype = barcode)


# function ----------------------------------------------------------------


# body --------------------------------------------------------------------

d <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell/all_hetero_af.cell.GSE226602_GSM7080029.qs")

d |>
  tidyr::pivot_longer(
    cols = -c(celltype, barcode),
    names_to = "variant",
    values_to = "af"
  ) |>
  tidyr::nest(
    .by = "variant",
    .key = "celltype_af"
  ) ->
d_nest

d_nest |>
  # head(100) |>
  dplyr::mutate(
    ks_test = parallel::mclapply(
      X = celltype_af,
      FUN = function(.x) {
        .x |>
          dplyr::filter(af > 0) |>
          dplyr::filter(celltype != "other") -> .xx
        tryCatch(
          {
            kruskal.test(
              af ~ celltype,
              data = .xx
            ) |>
              broom::tidy()
          },
          error = \(e){
            return(tibble::tibble(
              statistic = NA_real_,
              p.value = NA_real_,
              method = "Kruskal-Wallis test",
              parameter = NA_real_
            ))
          }
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(ks_test) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::arrange(p.value) ->
d_ks

d_ks$celltype_af[[1]] -> .x

.x |>
  dplyr::filter(celltype != "other") |>
  dplyr::filter(af > 0) |>
  ggplot(aes(
    x = af,
    y = celltype,
    fill = celltype
  )) +
  ggjoy::geom_joy()


FSA::dunnTest(
  af ~ celltype,
  data = .x,
  method = "bh"
) -> .dunn

.dunn$res |>
  dplyr::filter(P.adj < 0.05) |>
  dplyr::arrange(P.adj)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
