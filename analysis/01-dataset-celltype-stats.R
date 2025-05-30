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
sex_pred <- import("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir_sex.qs") |>
  dplyr::select(
    srrid,
    sex_pred = sex
  )

gse_dataset_metadata_full <- import(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.qs"
) |>
  dplyr::left_join(
    sex_pred,
    by = "srrid"
  ) |>
  dplyr::mutate(
    Gender = sex_pred
  )




# thegseid <- "GSE168453"
# body --------------------------------------------------------------------
gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.qs"
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
source("/home/liuc9/github/scMOCHA-data/stats/stats/00-colors.R")

gse_data_haplo_variant |>
  # dplyr::mutate(
  #   disease = dplyr::case_when(
  #     disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
  #     TRUE ~ "Other"
  #   )
  # ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      # levels = c(
      #   "Alzheimer's Disease",
      #   "COVID-19",
      #   "Healthy",
      #   "Unknown",
      #   "Other"
      # )
      levels = names(color_disease)
    )
  ) |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      # levels = c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE")
      levels = names(color_chemistry)
    )
  ) |>
  dplyr::arrange(disease, Chemistry) |>
  tidyr::unnest(cols = celltype_ratio) |>
  dplyr::mutate(
    celltype = factor(
      celltype,
      # levels = c(
      #   "B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T"
      # )
      levels = names(color_celltype)
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
  dplyr::arrange(disease, -b_ratio) |>
  dplyr::select(-data) ->
rank_srrid

# ggsci::pal_nejm()(5) |> color()


# disease_colors <- c(
#   "Alzheimer's Disease" = "#BC3C29FF",
#   "COVID-19" = "#0072B5FF",
#   "Healthy" = "#E18727FF",
#   "Unknown" = "grey50",
#   "Other" = "grey"
# )

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
  dplyr::ungroup() |>
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
  scale_fill_manual(
    name = "Disease",
    values = color_disease
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
  ) +
  scale_x_discrete(
    limits = rank_srrid$srrid,
  ) +
  theme(
    panel.background = element_blank(),
    # panel.background = element_rect(fill = "red", color = "black"),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 0, b = 0, unit = "cm")
  ) ->
p_tile_disease

# chemistry_colors <- c(
#   "SC5P-PE" = "#440154FF",
#   "SC5P-R2" = "#31688EFF",
#   "SC3Pv3" = "#35B779FF",
#   "SC3Pv2" = "#FDE725FF"
# )

rank_srrid |>
  dplyr::left_join(
    gse_dataset_metadata_full |>
      dplyr::select(srrid, Gender, Age_new, Chemistry),
    by = "srrid"
  ) |>
  dplyr::mutate(
    srrid = factor(srrid, levels = rank_srrid$srrid),
  ) |>
  ggplot(aes(
    x = srrid,
    y = 1
  )) +
  geom_tile(
    aes(
      fill = Chemistry
    )
  ) +
  scale_fill_manual(
    name = "Chemistry",
    values = color_chemistry
  ) +
  scale_x_discrete(
    limits = rank_srrid$srrid,
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 0, b = 0, unit = "cm")
  ) ->
p_tile_chemistry

# RColorBrewer::brewer.pal(8, "Set2") |> color()
# gender_colors <- c(
#   "Female" = "#FC8D62FF",
#   "Male" = "#66C2A5FF",
#   "Unknown" = "grey50"
# )

rank_srrid |>
  dplyr::mutate(
    srrid = factor(srrid, levels = rank_srrid$srrid),
  ) |>
  dplyr::left_join(
    gse_dataset_metadata_full |>
      dplyr::select(srrid, Gender, Age_new, Chemistry),
    by = "srrid"
  ) |>
  ggplot(aes(
    x = srrid,
    y = 1
  )) +
  geom_tile(
    aes(
      fill = Gender
    )
  ) +
  scale_fill_manual(
    name = "Sex",
    values = color_gender
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
  ) +
  scale_x_discrete(
    limits = rank_srrid$srrid,
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 0, b = 0, unit = "cm")
  ) ->
p_tile_gender

rank_srrid |>
  dplyr::mutate(
    srrid = factor(srrid, levels = rank_srrid$srrid),
  ) |>
  dplyr::left_join(
    gse_dataset_metadata_full |>
      dplyr::select(srrid, Gender, Age_new, Chemistry),
    by = "srrid"
  ) |>
  ggplot(aes(
    x = srrid,
    y = 1
  )) +
  geom_tile(
    aes(
      fill = Age_new
    )
  ) +
  scale_fill_gradient(
    name = "Age",
    low = "white",
    high = "red"
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
  ) +
  scale_x_discrete(
    limits = rank_srrid$srrid,
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 0, b = 0, unit = "cm")
  ) ->
p_tile_age

# celltype_colors <- c(
#   "B" = "#66C2A5FF",
#   "CD4 T" = "#FC8D62FF",
#   "CD8 T" = "#8DA0CBFF",
#   "DC" = "#E78AC3FF",
#   "Mono" = "#A6D854FF",
#   "NK" = "#FFD92FFF",
#   "other" = "#B3B3B3FF",
#   "other T" = "#E5C494FF"
# )

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
  # scale_fill_manual(
  #   name = "Cell Type",
  #   values = RColorBrewer::brewer.pal(8, "Set2")
  # ) +
  scale_fill_manual(
    name = "Cell Type",
    values = color_celltype
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
  ) +
  scale_x_discrete(
    limits = rank_srrid$srrid,
  ) +
  theme(
    plot.margin = margin(t = 0, b = 0, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black", ),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.y = element_text(
      color = "black",
      size = 18,
      face = "bold"
    ),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "Cell type ratio") ->
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
  # scale_fill_manual(
  #   name = "Cell Type",
  #   values = RColorBrewer::brewer.pal(8, "Set2")
  # ) +
  scale_fill_manual(
    name = "Cell Type",
    values = color_celltype
  ) +
  scale_y_continuous(
    expand = expansion(add = c(0.005, 0)),
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
    axis.title.y = element_text(
      color = "black",
      size = 18,
      face = "bold"
    ),
    axis.text.y = element_text(color = "black"),
  ) +
  labs(y = "# of cell") ->
p_celltype_count

{
  wrap_plots(
    p_celltype_ratio,
    p_celltype_count,
    p_tile_disease,
    p_tile_chemistry,
    p_tile_gender,
    p_tile_age,
    ncol = 1,
    heights = c(15, 15, 1, 1, 1, -1),
    guides = "collect"
  ) ->
  p_collect
  p_collect
}


ggsave(
  filename = file.path(outdir, "celltype_ratio.pdf" |> glue::glue()),
  plot = p_collect,
  width = 24,
  height = 12,
  dpi = 300
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
