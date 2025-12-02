#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-02 14:09:56
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)

gse_dataset_metadata_full <- import(
  "analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)
# load conn ---------------------------------------------------------------

# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

# body --------------------------------------------------------------------

gse_data |>
  dplyr::select(
    gseid,
    srrid,
    chemistry,
    anno,
    hetero,
    haplo_variant,
    haplo_violin,
    somatic_variant,
    celltype_ratio,
    clusteraf,
    bulkaf
  ) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(-gseid),
    by = c("srrid" = "srrid")
  ) |>
  dplyr::arrange(disease, Chemistry) -> gse_data_haplo_variant

gse_data_haplo_variant |>
  dplyr::mutate(
    heteroplasmic = purrr::map(
      somatic_variant,
      .f = \(.x) {
        # .x <- gse_data_haplo_variant$somatic_variant[[1]]
        # tibble::tibble(
        #   variant = .x$somatic
        # ) |>
        #   dplyr::mutate(
        #     pos = stringr::str_extract(variant, "\\d+") |> as.integer(),
        #   ) |>
        #   dplyr::filter(
        #     !pos %in% variants_tobe_excluded,
        #   ) -> .xx
        # .xx$variant -> heteroplasmic_variant
        # c(.x$high_af, .x$haplo) |> unique() -> homoplasmic_variant
        .x$heteroplasmic_variant <- .x$hete
        .x$homoplasmic_variant <- .x$homo

        .x
      }
    )
  ) |>
  dplyr::mutate(
    n_heteroplasmic = purrr::map(
      heteroplasmic,
      .f = \(.x) {
        tibble::tibble(
          n_heteroplasmic = length(.x$heteroplasmic_variant),
          n_homoplasmic = length(.x$homoplasmic_variant),
        )
      }
    )
  ) |>
  tidyr::unnest(cols = n_heteroplasmic) -> gse_data_variant_heteroplasmic


export(
  gse_data_variant_heteroplasmic,
  file = file.path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
    "gse_data_variant_heteroplasmic.qs"
  ),
  format = "qs",
)

#
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
