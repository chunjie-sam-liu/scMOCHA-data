#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-16 14:32:45
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

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"


dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
ks_test_dir <- file.path(dbdir, "all_hetero_af.cell.ks_test")
plotdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant"

all_heteroplasmic_af <- import(
  file.path(cleandatadir, "all_hetero_af.cluster.fst"),
) |>
  dplyr::select(-num_variants) |>
  dplyr::rename(celltype = barcode)


META <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.sex_pred.qs") |>
  dplyr::select(
    gseid, srrid, Age_new, Age_group,
    disease, Chemistry, sex_pred
  )


# function ----------------------------------------------------------------
fn_ks_test <- function(.gseid_srrid) {
  # .gseid_srrid <- "GSE226602_GSM7080017"
  .filename <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell/all_hetero_af.cell.{.gseid_srrid}.qs" |> glue::glue()
  d <- import(.filename)
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
        FUN = function(.m) {
          .m |>
            dplyr::filter(af > 0) |>
            dplyr::filter(celltype != "other") -> .mm
          tryCatch(
            {
              kruskal.test(
                af ~ celltype,
                data = .mm
              ) -> .fit
              .fit |>
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
        mc.cores = 20
      )
    ) |>
    tidyr::unnest(ks_test) |>
    dplyr::arrange(p.value) ->
  d_ks

  .filename_out <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/all_hetero_af.cell.ks_test/{.gseid_srrid}.ks_test.qs" |> glue::glue()
  export(
    d_ks,
    .filename_out,
  )
}


# body --------------------------------------------------------------------
all_heteroplasmic_af |>
  dplyr::select(gseid, srrid) |>
  dplyr::distinct() |>
  dplyr::mutate(
    gseid_srrid = paste(gseid, srrid, sep = "_")
  ) ->
gseid_srrid




# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
