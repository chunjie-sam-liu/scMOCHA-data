#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-16 18:18:19
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

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
expr <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.qs")

v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A.qs")



# body --------------------------------------------------------------------

expr |>
  dplyr::left_join(
    v_3173G_A,
    by = c("celltype")
  ) ->
expr_v_3173G_A

expr_v_3173G_A |>
  # head(20) |>
  dplyr::mutate(
    corr_results = parallel::mcmapply(
      function(.expr, .af) {
        # cor(x, y, method = "pearson")
        # .expr <- expr_v_3173G_A$expr[[1]]
        # .af <- expr_v_3173G_A$af[[1]]

        .expr |>
          dplyr::left_join(
            .af,
            by = c("gseid", "srrid")
          ) ->
        .expr_af
        tryCatch(
          {
            cor.test(~ af + expr, data = .expr_af, method = "pearson") |>
              broom::tidy() |>
              dplyr::select(
                corr = estimate,
                pval = p.value,
              )
          },
          error = function(e) {
            # message("Error in cor.test: ", e)
            return(NULL)
          }
        )
      },
      expr,
      af,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) |>
  dplyr::select(-expr, -af) |>
  tidyr::unnest(cols = corr_results) ->
expr_v_3173G_A_corr

export(
  expr_v_3173G_A_corr,
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A-corr.csv",
  format = "both"
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
