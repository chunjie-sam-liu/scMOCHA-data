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

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)

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
  dplyr::filter(gseid != "GSE220189") |>
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
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(cluster) |>
          dplyr::summarize(
            haplo_af_cluster = mean(af, na.rm = TRUE),
          ) ->
        .cluster_non_editing_cluster

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

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(cluster) |>
          dplyr::summarize(
            somatic_af_cluster = mean(af, na.rm = TRUE),
          ) ->
        .cluster_somatic_cluster

        tibble::tibble(
          haplo_af = mean(.variant_non_editing$mean_af, na.rm = TRUE),
          somatic_af = mean(.variant_somatic$mean_af, na.rm = TRUE),
          haplo_af_cell = list(.barcode_non_editing_cell),
          somatic_af_cell = list(.variant_somatic_cell),
          haplo_af_cluster = list(.cluster_non_editing_cluster),
          somatic_af_cluster = list(.cluster_somatic_cluster),
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
  geom_boxplot() +
  facet_wrap(
    ~chemistry,
    ncol = 3
  )


for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  # dplyr::filter(disease == "Healthy") |>
  ggplot(aes(
    x = Age_group,
    y = haplo_af,
  )) +
  geom_boxplot() +
  facet_wrap(
    ~chemistry,
    ncol = 3
  )

for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  # dplyr::filter(disease == "Healthy") |>
  ggplot(aes(
    x = Age_group,
    y = somatic_af,
  )) +
  # geom_boxplot() +
  geom_point() +
  facet_wrap(
    ~chemistry,
    ncol = 3
  )


outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"
celltypes <- c("B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T")


# ! somatic --------------------------------------------------------------------


for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  tidyr::unnest(cols = somatic_af_cluster) |>
  dplyr::mutate(
    cluster = factor(cluster, levels = celltypes)
  ) ->
for_hetero_af_forplot_cluster

for_hetero_af_forplot_cluster |>
  dplyr::group_by(Age_group) |>
  dplyr::summarise(
    mean_somatic_af = mean(somatic_af_cluster, na.rm = TRUE)
  ) ->
cluster_mean

for_hetero_af_forplot_cluster |>
  dplyr::left_join(
    cluster_mean,
    by = c("Age_group" = "Age_group")
  ) |>
  ggplot(aes(x = Age_group)) +
  ggh4x::facet_wrap2(
    ~cluster,
    ncol = 1,
    strip.position = "right",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = pcc$color
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    ),
    scales = "free_y",
  ) +
  geom_violin(
    aes(
      y = somatic_af_cluster,
      fill = mean_somatic_af
    ),
    alpha = 0.5,
    size = 1,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
  ) +
  ggnewscale::new_scale_fill() +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = Age_group,
      y = somatic_af_cluster,
      color = somatic_af_cluster
    ),
    size = 1,
    dodge.width = .75,
    alpha = .5,
    varwidth = TRUE
  ) +
  scale_color_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    limits = c(0, 1),
  ) +
  theme(
    plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "Mean heteroplasmy count"
  ) ->
p_age_somatic_af_cluster


ggsave(
  file.path(outdir, "p_age_somatic_af_cluster.pdf"),
  p_age_somatic_af_cluster,
  width = 12,
  height = 7
)



# ! haplo --------------------------------------------------------------------

for_hetero_af_forplot |>
  dplyr::filter(Age_group != "Unknown") |>
  dplyr::filter(!is.na(somatic_af)) |>
  dplyr::filter(!is.na(Age_new)) |>
  tidyr::unnest(cols = haplo_af_cluster) |>
  dplyr::mutate(
    cluster = factor(cluster, levels = celltypes)
  ) ->
for_hetero_af_forplot_cluster

for_hetero_af_forplot_cluster |>
  dplyr::group_by(Age_group) |>
  dplyr::summarise(
    mean_haplo_af = mean(haplo_af_cluster, na.rm = TRUE)
  ) ->
cluster_mean

for_hetero_af_forplot_cluster |>
  dplyr::left_join(
    cluster_mean,
    by = c("Age_group" = "Age_group")
  ) |>
  ggplot(aes(x = Age_group)) +
  ggh4x::facet_wrap2(
    ~cluster,
    ncol = 1,
    strip.position = "right",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = pcc$color
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      by_layer_y = FALSE,
    ),
    scales = "free_y",
  ) +
  geom_violin(
    aes(
      y = haplo_af_cluster,
      fill = mean_haplo_af
    ),
    alpha = 0.5,
    size = 1,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
  ) +
  ggnewscale::new_scale_fill() +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = Age_group,
      y = haplo_af_cluster,
      color = haplo_af_cluster
    ),
    size = 1,
    dodge.width = .75,
    alpha = .3,
    varwidth = TRUE
  ) +
  scale_color_gradient2(
    name = "AF",
    low = "white",
    mid = "red",
    high = "#3B0049",
    midpoint = 0.5,
  ) +
  scale_y_continuous(
    expand = c(0.01, 0),
    limits = c(0, 1),
  ) +
  theme(
    plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    # legend.position = c(0.8, 0.5),
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(
      size = 14,
      color = "black"
    ),
    legend.title = element_text(
      size = 16,
      colour = "black"
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      color = "black",
      face = "bold"
    )
  ) +
  labs(
    y = "Mean heteroplasmy count"
  ) ->
p_age_haplot_af_cluster


ggsave(
  file.path(outdir, "p_age_haplot_af_cluster.pdf"),
  p_age_haplot_af_cluster,
  width = 13,
  height = 7
)


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
