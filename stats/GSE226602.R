#!/usr/bin/env Rscript --vanilla
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: `r date()`
# @DESCRIPTION: filename

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
datadir <- "/home/liuc9/github/scMOCHA-data/data/GSE226602"
outdir <- file.path(datadir, "out")

# body --------------------------------------------------------------------

anno <- readr::read_rds(file.path(outdir, "GSE226602.scmocha.out.rds.gz"))
meta <- data.table::fread(file.path(outdir, "GSE226602.pheno.select.csv"))

anno |>
  dplyr::select(-c(srrdir, dir_exists)) |>
  dplyr::inner_join(meta, by = "srrid") ->
anno_meta

anno_meta |>
  dplyr::mutate(
    nmut = purrr::map_int(
      .x = anno,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_integer_)
        }
        nrow(.x)
      }
    )
  ) |>
  dplyr::mutate(
    nmut_somatic = purrr::map_int(
      .x = somatic_variant,
      .f = \(.x) {
        if (is.null(.x$somatic)) {
          return(NA_integer_)
        }
        length(.x$somatic)
      }
    )
  ) |>
  dplyr::mutate(
    haplogroup = purrr::map2(
      .x = anno,
      .y = srrid,
      .f = \(.x, .y) {
        message(.y)
        if (is.null(.x)) {
          return(
            tibble::tibble(
              Haplogroup = NA_character_,
              Verbose_haplogroup = NA_character_
            )
          )
        }
        .x |>
          dplyr::select(Haplogroup, Verbose_haplogroup) |>
          dplyr::filter(!is.na(Haplogroup)) |>
          dplyr::filter(Haplogroup != "") |>
          dplyr::distinct() |>
          dplyr::mutate_all(.funs = as.character) ->
        .xx

        if (nrow(.xx) == 0) {
          tibble::tibble(
            Haplogroup = NA_character_,
            Verbose_haplogroup = NA_character_
          )
        } else {
          .xx
        }
      }
    )
  ) |>
  tidyr::unnest(cols = haplogroup) |>
  tidyr::unnest(
    cols = cell_stats
  ) |>
  dplyr::mutate(
    ratio = round(`number of cells after filtering` / `estimated number of cells`, 2)
  ) ->
anno_meta_info

anno_meta_info |>
  dplyr::arrange(disease, age, gender) |>
  dplyr::select(
    Sample = srrid,
    Haplogroup = Haplogroup,
    Age = age,
    Gender = gender,
    Disease = disease,
    `# cells after filter` = `number of cells after filtering`,
    `# of variants` = nmut,
    `# of somatic variants` = nmut_somatic,
  ) ->
metadata_clean

# ggstats correlation plot ------------------------------------------------

ggstatsplot::ggscatterstats(
  data = metadata_clean |> dplyr::filter(Disease == "Healthy Control"),
  x = Age,
  y = `# of somatic variants`,
  title = "Correlation between Age and Number of Somatic Variants",
  xlab = "Age",
  ylab = "Number of Somatic Variants",
  ggtheme = ggplot2::theme_minimal()
) -> p
p


# ggsave(
#   filename = file.path(outdir, "correlation_age_somatic_variants.png"),
#   plot = p,
#   width = 8,
#   height = 6
# )
# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
