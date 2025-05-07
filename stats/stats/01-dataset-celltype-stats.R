#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-07 12:45:44
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

gse_dataset_metadata_full <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.rds"
)


pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)


# thegseid <- "GSE168453"
# body --------------------------------------------------------------------
gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.rds"
)


# body --------------------------------------------------------------------

gse_data |>
  dplyr::select(gseid, srrid, chemistry, anno, hetero, haplo_violin, somatic_variant) ->
for_hetero


gse_data |>
  dplyr::select(srrid, haplo_variant, hetero, celltype_ratio) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(srrid, Chemistry, disease),
    by = c("srrid" = "srrid")
  ) ->
gse_data_haplo_variant

# ! cell type ratio --------------------------------------------------------------------


gse_data_haplo_variant |>
  dplyr::mutate(
    disease = dplyr::case_when(
      disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
      TRUE ~ "Other"
    )
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Alzheimer's Disease",
        "COVID-19",
        "Healthy",
        "Unknown",
        "Other"
      )
    )
  ) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE")
    )
  ) |>
  dplyr::arrange(disease, Chemistry) |>
  tidyr::unnest(cols = celltype_ratio) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = c(
        "B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T"
      )
    )
  ) ->
for_celltype_ratio_plot

for_celltype_ratio_plot |>
  dplyr::group_by(srrid, disease) |>
  dplyr::select(srrid, disease, celltype, ratio) |>
  tidyr::nest() |>
  dplyr::mutate(
    b_ratio = purrr::map_dbl(
      .x = data,
      .f = \(.data) {
        .data |>
          dplyr::filter(celltype == "B") |>
          dplyr::pull(ratio) ->
        .b_ratio
        if (length(.b_ratio) == 0) {
          return(0)
        } else {
          return(.b_ratio)
        }
      }
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(disease, -b_ratio) ->
rank_srrid

rank_srrid |>
  dplyr::mutate(
    srrid = factor(srrid, levels = rank_srrid$srrid),
  ) |>
  dplyr::arrange(srrid) |>
  dplyr::group_by(disease) |>
  tibble::rowid_to_column() |>
  dplyr::mutate(
    mid_srrid = srrid[ceiling(dplyr::n() / 2)]
  ) |>
  ggplot(aes(
    x = srrid,
    y = 1
  )) +
  geom_tile(
    aes(
      fill = disease
    )
  ) +
  geom_text(
    aes(
      y = 1,
      label = ifelse(srrid == mid_srrid, as.character(disease), "")
    ),
  ) +
  ggsci::scale_fill_jco() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) ->
p_tile


for_celltype_ratio_plot |>
  dplyr::mutate(
    srrid = factor(
      srrid,
      levels = rank_srrid$srrid
    )
  ) |>
  ggplot(aes(
    x = srrid,
    y = ratio,
  )) +
  geom_col(
    aes(
      fill = celltype
    ),
    position = "stack"
  ) +
  scale_fill_manual(
    name = "Cell Type",
    values = RColorBrewer::brewer.pal(8, "Set2")
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0.005)),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "cell ratio") ->
p_celltype_ratio

for_celltype_ratio_plot |>
  dplyr::mutate(
    srrid = factor(
      srrid,
      levels = rank_srrid$srrid
    )
  ) |>
  ggplot(aes(
    x = srrid,
    y = n,
  )) +
  geom_col(
    aes(
      fill = celltype
    ),
    position = "stack"
  ) +
  scale_fill_manual(
    name = "Cell Type",
    values = RColorBrewer::brewer.pal(8, "Set2")
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0.005)),
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "# of cells") ->
p_celltype_count

ggsave(
  filename = file.path(outdir, "celltype_ratio.pdf" |> glue::glue()),
  plot = wrap_plots(
    p_celltype_ratio,
    p_celltype_count,
    p_tile,
    ncol = 1,
    heights = c(15, 15, 1),
    guides = "collect"
  ),
  width = 24,
  height = 12,
  dpi = 300
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
