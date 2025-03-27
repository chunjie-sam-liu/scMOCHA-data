#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-27 11:27:43
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
foundation_out <- file.path(basedir, "scfoundation/out")

gse_dataset_metadata_full <- readr::read_rds(
  file.path(foundation_out, "gse_dataset_metadata_full.rds")
)
gseids <- c(
  "GSE155673",
  "GSE157344",
  "GSE149689",
  "GSE171555",
  "GSE155223",
  "GSE163668",
  "GSE175524",
  "GSE206283",
  "GSE226598",
  "GSE261140",
  "GSE279945",
  "GSE214865",
  "GSE220189",
  "GSE233844",
  "GSE175499",
  "GSE149313",
  "GSE154386",
  "GSE159117",
  "GSE188632",
  "GSE166992",
  "GSE162117",
  "GSE226602",
  "GSE161354",
  "GSE235050",
  "GSE181279",
  # scfoundation2
  "GSE143353",
  "GSE148215",
  "GSE163314",
  "GSE163633",
  "GSE164690",
  "GSE167825",
  "GSE174125",
  "GSE184703",
  "GSE153421",
  "GSE147794",
  "GSE168453"
)


# body --------------------------------------------------------------------
tibble::tibble(
  gseid = gseids
) |>
  dplyr::mutate(
    anno = purrr::map(
      .x = gseid,
      .f = \(.gseid) {
        .anno <- readr::read_rds(
          file.path(basedir, .gseid, "out", glue::glue("{.gseid}.scmocha.out.rds.gz"))
        )
      }
    )
  ) ->
gse_data_loaded


gse_data_loaded |>
  tidyr::unnest(cols = anno) ->
gse_data

# body --------------------------------------------------------------------
gse_data
gse_data$haplo_violin[[1]]
gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_violin, somatic_variant) ->
for_hetero

# for_hetero$somatic_variant[[1]] -> .somatic_variant

for_hetero |>
  dplyr::mutate(
    mean_heteroplasmy_count = purrr::map2(
      .x = haplo_violin,
      .y = somatic_variant,
      .f = \(.haplo_violin, .somatic_variant) {
        .editing <- .somatic_variant$editing
        .somatic <- .somatic_variant$somatic
        .haplo_violin |>
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(variant) |>
          dplyr::summarize(
            mean_af = mean(af, na.rm = TRUE),
          ) ->
        .variant_non_editing

        .haplo_violin |>
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(barcode) |>
          dplyr::summarize(
            haplo_af_cell = mean(af, na.rm = TRUE),
          ) ->
        .barcode_non_editing_cell

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(variant) |>
          dplyr::summarize(
            mean_af = mean(af, na.rm = TRUE),
          ) ->
        .variant_somatic

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(barcode) |>
          dplyr::summarize(
            somatic_af_cell = mean(af, na.rm = TRUE),
          ) ->
        .variant_somatic_cell

        tibble::tibble(
          haplo_af = mean(.variant_non_editing$mean_af, na.rm = TRUE),
          somatic_af = mean(.variant_somatic$mean_af, na.rm = TRUE),
          haplo_af_cell = list(.barcode_non_editing_cell),
          somatic_af_cell = list(.variant_somatic_cell),
        )
      }
    )
  ) ->
for_hetero_af

for_hetero_af |>
  dplyr::select(gseid, srrid, chemistry, mean_heteroplasmy_count) |>
  tidyr::unnest(cols = mean_heteroplasmy_count) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(srrid, Age, Age_new, Age_group, disease),
    by = c("srrid" = "srrid")
  ) ->
for_hetero_af_forplot


for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(haplo_af)) |>
  ggplot(aes(
    x = Age_group,
    y = haplo_af,
  )) +
  geom_boxplot()

for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  # dplyr::filter(disease == "Healthy") |>
  ggplot(aes(
    x = Age_group,
    y = haplo_af,
  )) +
  geom_boxplot()

for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  # dplyr::filter(disease == "Healthy") |>
  ggplot(aes(
    x = Age_group,
    y = somatic_af,
  )) +
  geom_boxplot()


for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  tidyr::unnest(cols = somatic_af_cell) |>
  ggplot(aes(
    x = Age_group,
    y = mean_af,
  )) +
  geom_boxplot()


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
