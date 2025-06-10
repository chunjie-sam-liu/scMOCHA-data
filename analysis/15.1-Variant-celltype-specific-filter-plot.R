#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-06 10:58:26
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


# future: :plan(future: :multisession, workers = 10)


# load data ---------------------------------------------------------------
cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
dbdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/db"
ks_test_dir <- file.path(dbdir, "all_hetero_af.cell.ks_test")
plotdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-celltype-specific-variant"

# gseid_srrid_ks_load <- import(
#   file.path(
#     ks_test_dir,
#     "a_gseid_srrid_ks_load.qs"
#   )
# )

# gseid_srrid_ks_load |>
#   dplyr::filter(
#     p.value < 0.05,
#     statistic > 25
#   ) ->
# gseid_srrid_ks_load_p0.05_s55
# export(
#   gseid_srrid_ks_load_p0.05_s55,
#   file.path(
#     ks_test_dir,
#     "b_gseid_srrid_ks_load_p0.05_s25.qs"
#   )
# )

gseid_srrid_ks_load_p0.05_s55 <- import(
  file.path(
    ks_test_dir,
    "b_gseid_srrid_ks_load_p0.05_s25.qs"
  )
)

# function ----------------------------------------------------------------

fn_plot_joy <- function(.d, .variant = NULL) {
  # thevariant <- "7833T>C"
  if (is.null(.variant)) {
    .variant <- .d$variant[1]
  }

  .d |>
    dplyr::filter(af > 0) |>
    dplyr::mutate(
      celltype = gsub(
        "_",
        " ",
        celltype
      )
    ) |>
    dplyr::mutate(
      celltype = factor(
        celltype,
        levels = names(color_celltype) |> rev()
      )
    ) |>
    ggplot(aes(
      x = af,
      y = celltype,
      fill = celltype
    )) +
    ggjoy::geom_joy(
      # scale = 0.9,
      # alpha = 0.8,
      rel_min_height = 0.01,
      size = 0.1
    ) +
    scale_fill_manual(
      values = color_celltype,
      na.value = "grey50"
    ) +
    ggjoy::theme_joy() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        size = 16
      ),
    ) +
    labs(
      title = paste0(.variant),
      x = "Allele Frequency",
      y = "Cell Type"
    )
}
# body --------------------------------------------------------------------



# ? plot ks statistic--------------------------------------------------------------------

source("./analysis/00-colors.R")

gseid_srrid_ks_load |>
  dplyr::filter(p.value < 0.05) |>
  ggplot(
    aes(x = statistic)
  ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 100,
    fill = "grey50",
    color = "black",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 55,
    linetype = "dashed",
    color = "red"
  ) +
  geom_text(
    data = data.frame(
      x = 55,
      y = 0.02,
      label = "KS statistic = 55",
      vjust = -1
    ),
    aes(
      x = 55,
      y = 0.02,
      label = "KS statistic = 55",
      vjust = -1
    ),
    color = "red",
    size = 4
  ) +
  theme_bw() +
  labs(
    x = "KS statistic",
    y = "Density",
    title = "Distribution of KS statistic for all variants"
  ) ->
plot_ks_statistic
ggsave(
  file.path(
    plotdir,
    "ks_statistic_distribution.pdf"
  ),
  plot = plot_ks_statistic,
  width = 8,
  height = 6
)



# ? find examples --------------------------------------------------------------------

gseid_srrid_ks_load |>
  dplyr::filter(
    p.value < 0.05,
    statistic > 20
  ) |>
  dplyr::group_by(
    variant
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_statistic = mean(statistic, na.rm = TRUE),
    mean_p_value = mean(p.value, na.rm = TRUE),
  ) |>
  dplyr::mutate(
    mean_log10p = -log10(mean_p_value),
  ) ->
variant_count_statistic

variant_count_statistic$mean_log10p |> summary()
variant_count_statistic$mean_statistic |> summary()

variant_count_statistic |>
  ggplot(aes(
    x = mean_statistic,
    y = mean_log10p,
    size = n
  )) +
  ggpointdensity::geom_pointdensity(
    adjust = 0.01,
    # show.legend = FALSE,
  ) +
  # geom_point() +
  viridis::scale_color_viridis() +
  scale_x_continuous(
    limits = c(20, 170),
    labels = scales::number,
    breaks = seq(20, 170, by = 10),
    expand = expansion(add = c(0.01, 0))
  ) +
  scale_y_continuous(
    limits = c(2, 15),
    labels = scales::number,
    breaks = seq(2, 15, by = 1),
    expand = expansion(add = c(0.01, 0))
  ) +
  ggrepel::geom_text_repel(
    data = variant_count_statistic |>
      dplyr::filter(
        variant %in% c(
          "7833T>C",
          # sort by mean_statistic
          "3030A>G",
          "7430A>C",
          "4175G>A",
          "3727T>C",
          "7418C>A",
          "6409T>C",
          "6669C>G",
          "7583T>G",
          "7582C>G",
          "929A>C",
          # sort by mean_log10p
          "4886C>T",
          "1082A>G",
          # "15213T>C",
          "3173G>A",
          "6669C>G",
          "13956A>G",
          "14063T>C",
          "8849T>C",
          "8285C>A",
          "4794G>A",
          "11502T>C",
          # sort by n
          "4175G>A",
          "2289G>T",
          "10645T>G",
          "3584A>C",
          "3030A>G",
          "3520A>C",
          "2193T>A",
          "9076A>C",
          "8072T>G",
          "3577A>C"
        )
      ),
    aes(label = variant),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE,
    nudge_x = 5,
    nudge_y = 2
  ) +
  labs(
    x = "Mean KS Statistic",
    y = "-log10(Mean P-value)",
    title = "Variant Count vs. Mean KS Statistic and P-value"
  ) +
  theme_bw() +
  guides(
    size = guide_legend(title = "Variant Count"),
    color = guide_colorbar(title = "Density")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
  ) ->
plot_variant_count_statistic
ggsave(
  file.path(
    plotdir,
    "ks_variant_count_statistic.pdf"
  ),
  plot = plot_variant_count_statistic,
  width = 11,
  height = 6
)




# ? example --------------------------------------------------------------------
variant_list <- c(
  "7833T>C",
  # sort by mean_statistic
  "3030A>G",
  "7430A>C",
  "4175G>A",
  "3727T>C",
  "7418C>A",
  "6409T>C",
  "6669C>G",
  "7583T>G",
  "7582C>G",
  "929A>C",
  # sort by mean_log10p
  # "4886C>T",
  "1082A>G",
  # "15213T>C",
  "3173G>A",
  "6669C>G",
  "13956A>G",
  "14063T>C",
  "8849T>C",
  "8285C>A",
  "4794G>A",
  "11502T>C",
  # sort by n
  "4175G>A",
  "2289G>T",
  "10645T>G",
  "3584A>C",
  "3030A>G",
  "3520A>C",
  "2193T>A",
  "9076A>C",
  "8072T>G",
  "3577A>C"
)
thevariant <- "7833T>C"
thevariant <- "3727T>C"

variant_count_statistic |>
  dplyr::arrange(-mean_statistic) |>
  dplyr::slice(1:100) |>
  dplyr::pull(variant) ->
variant_list


length(variant_list)
variant_list |>
  purrr::map(function(thevariant) {
    gseid_srrid_ks_load |>
      dplyr::filter(variant == thevariant) |>
      dplyr::filter(
        p.value < 0.05,
        # statistic > 100
      ) |>
      dplyr::slice(1) |>
      tidyr::unnest(cols = celltype_af) |>
      fn_plot_joy()
  }) ->
plot_variants_list

plot_variants_list |>
  wrap_plots(ncol = 10) +
  plot_layout(
    guides = "collect",
  ) ->
plot_variants_joy
ggsave(
  plot = plot_variants_joy,
  filename = file.path(
    plotdir,
    "joy_variants100.pdf"
  ),
  width = 50,
  height = 50,
  limitsize = FALSE
)



# ? single variant joy --------------------------------------------------------------------

thevariant <- "3173G>A"
thevariant <- "7833T>C"
thevariant <- "3727T>C"
gseid_srrid_ks_load |>
  dplyr::filter(variant == thevariant) |>
  dplyr::filter(
    p.value < 0.05,
    # statistic > 100
  ) |>
  dplyr::mutate(
    p = purrr::map2(
      .x = gseid_srrid,
      .y = celltype_af,
      .f = function(.x, .y) {
        fn_plot_joy(
          .d = .y,
          .variant = .x
        )
      }
    )
  ) ->
plot_thevariant_sample_list


plot_thevariant_sample_list |>
  dplyr::slice(1:100) |>
  dplyr::pull(p) |>
  wrap_plots(ncol = 10) +
  plot_layout(
    guides = "collect",
  ) ->
plot_thevariant_sample_joy
ggsave(
  plot = plot_thevariant_sample_joy,
  filename = file.path(
    plotdir,
    "{thevariant}_joy_sample.pdf" |> glue::glue()
  ),
  width = 50,
  height = 50,
  limitsize = FALSE
)




# ? GSE226602_GSM7080017 azimuth.rda --------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
GSE226602_GSM7080017 <- import("/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE226602/final/GSM7080017/sc_azimuth.rds.gz")
GSE226602_GSM7080017$sc_azimuth@meta.data |> head()

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
