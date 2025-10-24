#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-02 10:55:16
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
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz"
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-basic"

gse_dataset_metadata_full <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)


gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data_fisher.qs"
)

gse_data_variant_heteroplasmic <- import(file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/",
  "gse_data_variant_heteroplasmic_fisher.qs"
)) |>
  dplyr::select(gseid, srrid, n_heteroplasmic, n_homoplasmic)


gse_dataset_metadata_full |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(
    samples = dplyr::n(),
    `Avg. Age` = mean(Age_new, na.rm = T),
    `Avg. # of cells` = mean(`# cells after filter`, na.rm = T),
    `Avg. median genes/cell` = mean(`Median genes/cell`, na.rm = T),
    `Avg. median UMI/cell` = mean(`Median UMI/cell`, na.rm = T),
    Chemistry = paste0(sort(unique(Chemistry)), collapse = ", "),
    Disease = paste0(sort(unique(disease)), collapse = ", "),
    Source = paste0(sort(unique(Source)), collapse = ", "),
    Publication = paste0(sort(unique(Publication)), collapse = ", "),
  ) -> gse_dataset_metadata_full_sel

gse_data |>
  dplyr::select(
    gseid,
    srrid,
    metrics,
    depth_read,
    depth,
    somatic_variant = somatic_variant_fisher
  ) |>
  dplyr::left_join(
    gse_data_variant_heteroplasmic,
    by = c("gseid", "srrid")
  ) |>
  dplyr::mutate(
    total_reads = purrr::map_dbl(
      .x = metrics,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        .x$`Number of Reads`
      }
    ),
    depth_read_mean = purrr::map_dbl(
      .x = depth_read,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        mean(.x$depth, na.rm = T)
      }
    ),
    depth_mean = purrr::map_dbl(
      .x = depth,
      .f = \(.x) {
        if (is.null(.x)) {
          return(NA_real_)
        }
        mean(.x$depth, na.rm = T)
      }
    ),
  ) |>
  dplyr::mutate(
    nmut_variant = purrr::map2(
      # .x = somatic_variant,
      .x = n_heteroplasmic,
      .y = n_homoplasmic,
      .f = \(.x, .y) {
        # .x <- a$somatic_variant[[1]]
        # union(.x$haplo, .x$somatic) |>
        #   unique() |>
        #   length() -> .nmut
        # .x |>
        #   purrr::map_int(length) |>
        #   tibble::enframe() |>
        #   tidyr::spread(key = name, value = value) -> .xx
        # names(.xx) <- glue::glue("nmut_{names(.xx)}")
        # .xx |>
        #   tibble::add_column(
        #     nmut = .nmut,
        #     .before = 1
        #   )
        tibble::tibble(
          nmut = sum(.x, .y, na.rm = T),
          nmut_somatic = .x,
          nmut_haplo = .y,
        )
      }
    )
  ) |>
  tidyr::unnest(cols = nmut_variant) -> gse_data_read

# body --------------------------------------------------------------------
source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()

gse_data_read |>
  dplyr::select(
    gseid,
    nmut,
    nmut_somatic,
    total_reads,
    depth_read_mean,
    depth_mean
  ) |>
  dplyr::group_by(gseid) |>
  dplyr::summarise(
    `Avg. somatic mutation` = mean(nmut_somatic, na.rm = T),
    `Avg. mutation` = mean(nmut, na.rm = T),
    `Avg. mapped reads` = mean(depth_read_mean, na.rm = T),
    `Avg. total reads` = mean(total_reads, na.rm = T),
    # `Avg. call depth` = mean(depth_mean, na.rm = T),
  ) |>
  dplyr::left_join(
    gse_dataset_metadata_full_sel,
    by = "gseid"
  ) |>
  dplyr::relocate(
    gseid,
    samples,
    Chemistry,
    .before = 1
  ) |>
  dplyr::rename(
    `GSE ID` = gseid,
    Samples = samples,
  ) -> gses_meta_read_all

gses_meta_read_all |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = names(color_chemistry)),
  ) |>
  dplyr::arrange(
    Chemistry,
    -`Avg. somatic mutation`
  ) -> df

# chem_colors <- viridis::viridis_pal(option = "D")(4) |>
#   prismatic::color()

# df |> dplyr::glimpse()

the_header <- data.frame(
  col_keys = c(
    "GSE ID",
    "Samples",
    "Chemistry",
    "Avg. somatic mutation",
    "Avg. mutation",
    "Avg. mapped reads",
    "Avg. total reads",
    "Avg. Age",
    "Avg. # of cells",
    "Avg. median genes/cell",
    "Avg. median UMI/cell",
    "Disease",
    "Source",
    "Publication"
  ),
  line2 = c(
    "GSE ID",
    "Samples",
    "Chemistry",
    "# of mutations",
    "# of mutations",
    "# of reads",
    "# of reads",
    "Avg. Age",
    "Avg. # of cells",
    "Avg. median genes/cell",
    "Avg. median UMI/cell",
    "Disease",
    "Source",
    "Publication"
  ),
  line3 = c(
    "GSE ID",
    "Samples",
    "Chemistry",
    "Somatic",
    "Total",
    "mtDNA",
    "Total",
    "Avg. Age",
    "Avg. # of cells",
    "Avg. median genes/cell",
    "Avg. median UMI/cell",
    "Disease",
    "Source",
    "Publication"
  )
)

chem_uniq <- df$Chemistry |>
  unique() |>
  as.character()
last_indices <- sapply(chem_uniq, function(x) {
  max(which(df$Chemistry == x))
})
# RColorBrewer::brewer.pal(5, "Set2") |> prismatic::color()

library(flextable)

the_header |>
  tibble::rowid_to_column() |>
  dplyr::group_by(line2) |>
  dplyr::filter(!dplyr::n() > 1) -> the_header_idx

flextable::flextable(df) |>
  flextable::set_header_df(
    the_header,
    key = "col_keys"
  ) |>
  flextable::merge_v(part = "header", j = the_header_idx$rowid) |>
  flextable::merge_h(part = "header", i = c(1, 2)) |>
  flextable::theme_booktabs(bold_header = TRUE) |>
  flextable::bg(
    i = ~ Chemistry == names(color_chemistry)[1],
    j = c("Chemistry"),
    bg = color_chemistry[1]
  ) |>
  flextable::bg(
    i = ~ Chemistry == names(color_chemistry)[2],
    j = c("Chemistry"),
    bg = color_chemistry[2]
  ) |>
  flextable::bg(
    i = ~ Chemistry == names(color_chemistry)[3],
    j = c("Chemistry"),
    bg = color_chemistry[3]
  ) |>
  flextable::bg(
    i = ~ Chemistry == names(color_chemistry)[4],
    j = c("Chemistry"),
    bg = color_chemistry[4]
  ) |>
  flextable::bg(
    bg = scales::col_numeric(
      palette = c("transparent", "#FFFFB3FF"),
      domain = c(min(df$`Avg. somatic mutation`), max(df$`Avg. mutation`))
    ),
    j = c(4, 5),
    part = "body"
  ) |>
  flextable::bg(
    bg = scales::col_numeric(
      palette = c("transparent", "#FB8072FF"),
      domain = c(min(df$`Avg. mapped reads`), max(df$`Avg. mapped reads`))
    ),
    j = c(6),
    part = "body"
  ) |>
  flextable::bg(
    bg = scales::col_numeric(
      palette = c("transparent", "#FB8072FF"),
      domain = c(min(df$`Avg. total reads`), max(df$`Avg. total reads`))
    ),
    j = c(7),
    part = "body"
  ) |>
  flextable::color(
    j = c("Chemistry"),
    color = "white"
  ) |>
  flextable::bold(
    j = c("GSE ID", "Chemistry")
  ) |>
  flextable::italic(
    j = c("Publication")
  ) |>
  flextable::vline(
    j = c(
      "Samples",
      "Chemistry",
      "Avg. mutation",
      "Avg. total reads",
      "Avg. median UMI/cell"
    ),
    border = flextable::fp_border_default()
  ) |>
  flextable::hline(
    i = last_indices,
    border = flextable::fp_border_default()
  ) |>
  flextable::colformat_double(
    j = c("Avg. somatic mutation", "Avg. mutation", "Avg. Age")
  ) |>
  flextable::colformat_num(
    j = c(
      "Avg. total reads",
      "Avg. mapped reads",
      "Avg. # of cells",
      "Avg. median genes/cell",
      "Avg. median UMI/cell"
    ),
  ) |>
  flextable::align(align = "center", part = "all") |>
  flextable::valign(valign = "center", part = "header") |>
  flextable::width(
    j = c(6, 7, 12),
    width = 1.2
  ) |>
  flextable::width(
    j = c(14),
    width = 2
  ) -> ft
ft
flextable::save_as_image(
  ft,
  path = file.path(outdir, "gses_meta_read_fisher.svg"),
  width = 20,
  height = 7
)

flextable::save_as_pptx(
  ft,
  path = file.path(outdir, "gses_meta_read_fisher.pptx")
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
