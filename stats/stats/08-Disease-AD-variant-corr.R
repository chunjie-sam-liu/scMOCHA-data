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




# body --------------------------------------------------------------------



# !  v_3173G_A--------------------------------------------------------------------



v_3173G_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3173G>A.qs")
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


# expr_v_3173G_A_corr |>
#   dplyr::filter(celltype == "Mono") |>
#   dplyr::filter(pval < 0.05, abs(corr) > 0.5) |>
#   dplyr::arrange(desc(abs(corr)))



# ! v_1397T_A --------------------------------------------------------------------

v_1397T_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-1397T>A.qs")
expr |>
  dplyr::left_join(
    v_1397T_A,
    by = c("celltype")
  ) ->
expr_v_1397T_A


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
expr_v_1397T_A_corr

export(
  expr_v_1397T_A_corr,
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-1397T>A-corr.csv",
  format = "both"
)


# ! v_1670A_G --------------------------------------------------------------------

v_1670A_G <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-1670A>G.qs")
expr |>
  dplyr::left_join(
    v_1670A_G,
    by = c("celltype")
  ) ->
expr_v_1670A_G

expr_v_1670A_G |>
  dplyr::mutate(
    corr_results = parallel::mcmapply(
      function(.expr, .af) {
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
expr_v_1670A_G_corr

export(
  expr_v_1670A_G_corr,
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-1670A>G-corr.csv",
  format = "both"
)


# ! v_3176A_T --------------------------------------------------------------------

v_3176A_T <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3176A>T.qs")
expr |>
  dplyr::left_join(
    v_3176A_T,
    by = c("celltype")
  ) ->
expr_v_3176A_T

expr_v_3176A_T |>
  dplyr::mutate(
    corr_results = parallel::mcmapply(
      function(.expr, .af) {
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
expr_v_3176A_T_corr

export(
  expr_v_3176A_T_corr,
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3176A>T-corr.csv",
  format = "both"
)


# ! v_3178T_A --------------------------------------------------------------------

v_3178T_A <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3178T>A.qs")
expr |>
  dplyr::left_join(
    v_3178T_A,
    by = c("celltype")
  ) |>
  # Filter rows where af is not NULL to ensure proper size matching
  dplyr::filter(!purrr::map_lgl(af, is.null)) ->
expr_v_3178T_A

expr_v_3178T_A |>
  dplyr::mutate(
    corr_results = purrr::map2(
      expr, af,
      function(.expr, .af) {
        .expr |>
          dplyr::left_join(
            .af,
            by = c("gseid", "srrid")
          ) ->
        .expr_af

        # Only proceed if we have sufficient data after joining
        if (nrow(.expr_af) < 3 || all(is.na(.expr_af$af)) || all(is.na(.expr_af$expr))) {
          return(tibble::tibble(corr = NA_real_, pval = NA_real_))
        }

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
            return(tibble::tibble(corr = NA_real_, pval = NA_real_))
          }
        )
      }
    )
  ) |>
  dplyr::select(-expr, -af) |>
  tidyr::unnest(cols = corr_results) ->
expr_v_3178T_A_corr

export(
  expr_v_3178T_A_corr,
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/ad/ad-celltype-variant-af-3178T>A-corr.csv",
  format = "both"
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
