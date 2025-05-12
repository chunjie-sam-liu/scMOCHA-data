#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-09 16:26:12
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
log_threshold(TRACE)
log_layout(layout_glue_colors)

# future: :plan(future: :multisession, workers = 10)

# function ----------------------------------------------------------------


# load data ---------------------------------------------------------------
basedir <- "/home/liuc9/github/scMOCHA-data/data"
outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

gse_dataset_metadata_full <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.rds"
)

gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.rds"
)

all_variant <- readr::read_rds("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/all_variant.rds") |>
  dplyr::select(variant, issomatic)

all_heteroplasmic_af <- data.table::fread(
  file.path(cleandatadir, "all_hetero_af.bulk.csv"),
  header = TRUE,
  sep = ",",
) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, barcode, num_variants),
    names_to = "variant",
    values_to = "af"
  )

# body --------------------------------------------------------------------


gse_dataset_metadata_full |>
  dplyr::filter(
    disease %in% c("Healthy", "COVID-19")
  ) |>
  dplyr::select(
    gseid, srrid, Chemistry, disease
  ) ->
admeta

admeta |>
  dplyr::count(
    disease, Chemistry
  )

admeta |>
  dplyr::filter(
    Chemistry != "SC5P-PE"
  ) ->
admeta_sc5p


admeta_sc5p |>
  dplyr::left_join(
    gse_data,
    by = c("gseid", "srrid")
  ) ->
admeta_sc5p_variant

admeta_sc5p_variant |>
  dplyr::mutate(
    variant_type = purrr::map(
      .x = anno,
      .f = function(.x) {
        .x |>
          dplyr::mutate(
            variant = glue::glue("{Position}{Ref}>{Alt}")
          ) |>
          dplyr::select(variant, ntchange) |>
          dplyr::left_join(
            all_variant,
            by = "variant"
          )
      }
    )
  ) |>
  dplyr::select(
    gseid, srrid, Chemistry, disease, variant_type
  ) ->
admeta_sc5p_variant_type



# ! variant venn --------------------------------------------------------------------



admeta_sc5p_variant_type |>
  tidyr::unnest(cols = variant_type) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) |>
  dplyr::select(
    srrid, disease, variant
  ) |>
  dplyr::group_by(
    variant
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    m = purrr::map(
      .x = data,
      .f = function(.x) {
        .x |>
          dplyr::group_by(disease) |>
          dplyr::count() |>
          dplyr::ungroup() |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n
          )
      }
    )
  ) |>
  tidyr::unnest(cols = m) ->
admeta_sc5p_variant_type_count

admeta_sc5p_variant_type_count |>
  dplyr::select(-data) |>
  dplyr::arrange(
    Healthy
  ) |>
  dplyr::filter(
    variant == "3173G>A"
  )



# ! compare variant between AD and healthy --------------------------------------------------------------------

admeta_sc5p_variant_type |>
  tidyr::unnest(cols = variant_type) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) |>
  dplyr::select(
    srrid, disease, variant
  ) |>
  dplyr::left_join(
    all_heteroplasmic_af,
    by = c("srrid", "variant")
  ) ->
admeta_sc5p_variant_type_af

admeta_sc5p_variant_type_af |>
  dplyr::group_by(
    variant
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    m = purrr::map(
      .x = data,
      .f = \(.x) {
        # .x <- a$data[[1]]
        tryCatch(
          expr = {
            t.test(
              af ~ disease,
              data = .x,
              var.equal = TRUE
            ) |>
              broom::tidy() |>
              dplyr::select(
                estimate, estimate1, estimate2, p.value, conf.low, conf.high
              )
          },
          error = function(e) {
            message("Error: ", conditionMessage(e))
            return(NULL)
          }
        )
      }
    )
  ) ->
admeta_sc5p_variant_type_af_ttest

admeta_sc5p_variant_type_af_ttest |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = m) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    est = abs(estimate),
  ) |>
  dplyr::mutate(
    rank = plog10p * est,
  ) |>
  dplyr::arrange(
    desc(rank)
  ) ->
admeta_sc5p_variant_type_af_ttest_rank


admeta_sc5p_variant_type_af |>
  dplyr::filter(
    variant == admeta_sc5p_variant_type_af_ttest_rank$variant[1]
  ) |>
  ggstatsplot::ggbetweenstats(
    data = _,
    x = disease,
    y = af,
    pairwise.display = "p-value",
    pairwise.comparisons = TRUE,
    p.adjust.method = "fdr",
    p.adjust.display = TRUE,
    p.value.label = "p.adj",
    p.value.label.size = 3.5,
    p.value.label.color = "black",
    p.value.label.position = c(0.5, 0.95),
    p.value.label.nudge_y = 0.05,
    ggplot.component = list(
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    )
  )

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
