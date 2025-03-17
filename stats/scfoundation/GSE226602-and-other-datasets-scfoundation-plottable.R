#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-20 17:24:44
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

datadir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation/out"
filename_ <- "gses_meta_read.xlsx"
gses_meta_read <- readxl::read_xlsx(
  file.path(
    datadir,
    filename_
  )
)

datadir <- "/home/liuc9/github/scMOCHA-data/data/out_new_ting"
gses_meta_read_ <- readxl::read_xlsx(
  file.path(
    datadir,
    filename_
  )
)

datadir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC/out"
gses_meta_read__ <- readxl::read_xlsx(
  file.path(
    datadir,
    filename_
  )
)

# body --------------------------------------------------------------------
chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()

gses_meta_read |>
  dplyr::bind_rows(gses_meta_read_) |>
  dplyr::bind_rows(gses_meta_read__) ->
gses_meta_read_all

gses_meta_read_all$gseid |>
  sort() |>
  unique() -> gseids

gses_meta_read_all |>
  dplyr::filter(gseid %in% gseids) |>
  dplyr::distinct() |>
  dplyr::select(
    `GSE ID` = gseid,
    samples,
    Chemistry,
    `Avg. somatic mutation`,
    `Avg. mutation`,
    `Avg. mapped reads`,
    `Avg. total reads`,
    Disease,
    Source,
    Publication
  ) |>
  dplyr::mutate(
    Chemistry = factor(Chemistry, levels = chem_levels)
  ) |>
  dplyr::arrange(
    Chemistry, -`Avg. somatic mutation`
  ) ->
df


library(flextable)

ggsci::pal_aaas()(3) |> prismatic::color()
chem_colors <- viridis::viridis_pal(option = "D")(4) |>
  prismatic::color()

the_header <- data.frame(
  col_keys = c(
    "GSE ID",
    "samples",
    "Chemistry",
    "Avg. somatic mutation",
    "Avg. mutation",
    "Avg. mapped reads",
    "Avg. total reads",
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
df$samples |> sum()
flextable::flextable(df) |>
  flextable::set_header_df(
    the_header,
    key = "col_keys"
  ) |>
  flextable::merge_v(part = "header", j = c(1, 2, 3, 8, 9, 10)) |>
  flextable::merge_h(part = "header", i = c(1, 2)) |>
  theme_booktabs(bold_header = TRUE) |>
  flextable::bg(
    i = ~ Chemistry == chem_levels[1],
    j = c("Chemistry"),
    bg = chem_colors[1]
  ) |>
  flextable::bg(
    i = ~ Chemistry == chem_levels[2],
    j = c("Chemistry"),
    bg = chem_colors[2]
  ) |>
  flextable::bg(
    i = ~ Chemistry == chem_levels[3],
    j = c("Chemistry"),
    bg = chem_colors[3]
  ) |>
  flextable::bg(
    i = ~ Chemistry == chem_levels[4],
    j = c("Chemistry"),
    bg = chem_colors[4]
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
    j = c("Chemistry", "Avg. mutation", "Avg. total reads"),
    border = fp_border_default()
  ) |>
  flextable::hline(
    i = last_indices,
    border = fp_border_default()
  ) |>
  colformat_double(
    j = c("Avg. somatic mutation", "Avg. mutation")
  ) |>
  colformat_num(
    j = c("Avg. total reads", "Avg. mapped reads")
  ) |>
  align(align = "center", part = "all") |>
  valign(valign = "center", part = "header") |>
  flextable::width(
    j = c(7, 8),
    width = 1.2
  ) |>
  flextable::width(
    j = c(10),
    width = 2
  ) ->
ft
ft

datadir <- "/home/liuc9/github/scMOCHA-data/data/scfoundation/out"
flextable::save_as_image(
  ft,
  path = file.path(datadir, "gses_meta_read.svg"),
  width = 20,
  height = 7
)

flextable::save_as_pptx(
  ft,
  path = file.path(datadir, "gses_meta_read.pptx")
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
