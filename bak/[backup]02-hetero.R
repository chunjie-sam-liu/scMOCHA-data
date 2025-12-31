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
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz"

gse_dataset_metadata_full <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.rds"
)


pcc <- readr::read_tsv(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)


# thegseid <- "GSE168453"
# body --------------------------------------------------------------------
gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.rds"
)

# body --------------------------------------------------------------------

gse_data |>
  dplyr::select(
    gseid,
    srrid,
    chemistry,
    anno,
    hetero,
    haplo_violin,
    somatic_variant
  ) -> for_hetero


gse_data |>
  dplyr::select(srrid, haplo_variant, hetero, celltype_ratio) |>
  dplyr::left_join(
    gse_dataset_metadata_full |> dplyr::select(srrid, Chemistry, disease),
    by = c("srrid" = "srrid")
  ) -> gse_data_haplo_variant


# gse_data_haplo_variant |>
#   dplyr::select(haplo_variant) |>
#   tidyr::unnest(cols = haplo_variant) |>
#   dplyr::select(
#     Position, Ref, Alt, Locus, Haplogroup, Verbose_haplogroup, Disease, variant, color
#   ) |>
#   dplyr::distinct() ->
# all_variants

# all_variants |>
#   dplyr::filter(color == "white") ->
# haplo_variants
# all_variants |>
#   dplyr::filter(color == "black") ->
# somatic_variants

# ! no longer use --------------------------------------------------------------------

\() {
  gse_data_haplo_variant |>
    tidyr::unnest(cols = hetero) |>
    dplyr::filter(variant %in% somatic_variants$variant) |>
    dplyr::filter(disease %in% c("Alzheimer's Disease", "Healthy")) |>
    dplyr::group_by(
      srrid,
      disease,
      variant
    ) |>
    dplyr::summarise(af = mean(af, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::arrange(disease) -> for_hetero_af_tile

  for_hetero_af_tile |>
    dplyr::select(variant) |>
    dplyr::mutate(
      pos = gsub(variant, pattern = "[ACGT]|>", replacement = "")
    ) |>
    dplyr::mutate(pos = as.numeric(pos)) |>
    dplyr::arrange(-pos) |>
    dplyr::distinct() -> rank_variant

  for_hetero_af_tile |>
    dplyr::group_by(disease, srrid) |>
    dplyr::count() |>
    dplyr::ungroup() |>
    dplyr::arrange(disease, -n) -> rank_srrid

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
    ggsci::scale_fill_npg() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) -> p_tile

  for_hetero_af_tile |>
    dplyr::mutate(
      srrid = factor(
        srrid,
        levels = rank_srrid$srrid
      ),
      variant = factor(
        variant,
        levels = rank_variant$variant
      )
    ) |>
    ggplot(aes(
      x = srrid,
      y = variant,
      fill = disease
    )) +
    geom_tile(
      show.legend = F
    ) +
    ggsci::scale_fill_npg() +
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
    labs(
      y = "Variant",
    ) -> p_variants_tile

  ggsave(
    filename = file.path(outdir, "disease_varaint.pdf" |> glue::glue()),
    plot = wrap_plots(
      p_variants_tile,
      p_tile,
      ncol = 1,
      heights = c(15, 1),
      guides = "collect"
    ),
    width = 14,
    height = 7,
    dpi = 300
  )
}


# ! mean heteroplasmic variant count --------------------------------------------------------------------

all_variant <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.rds"
)

# for_hetero$somatic_variant[[1]] -> .somatic_variant
for_hetero |>
  dplyr::mutate(
    mean_heteroplasmy_count = parallel::mcmapply(
      FUN = \(.haplo_violin, .somatic_variant, i, total) {
        cli::cli_alert_info("Processing {i}/{total} haplo_violin...")
        # for_hetero$somatic_variant[[1]] -> .somatic_variant
        # for_hetero$haplo_violin[[1]] -> .haplo_violin

        .editing <- all_variant |>
          dplyr::filter(issomatic %in% c("multiple", "other")) |>
          dplyr::pull(variant)
        .somatic <- all_variant |>
          dplyr::filter(issomatic == "heteroplasmic") |>
          dplyr::pull(variant)

        .haplo_violin |>
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(variant) |>
          dplyr::summarize(
            mean_af = mean(af, na.rm = TRUE),
          ) -> .variant_non_editing

        .haplo_violin |>
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(barcode) |>
          dplyr::summarize(
            haplo_af_cell = mean(af, na.rm = TRUE),
          ) -> .barcode_non_editing_cell

        .haplo_violin |>
          dplyr::filter(!variant %in% .editing) |>
          dplyr::group_by(cluster) |>
          dplyr::summarize(
            haplo_af_cluster = mean(af, na.rm = TRUE),
          ) -> .cluster_non_editing_cluster

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(variant) |>
          dplyr::summarize(
            mean_af = mean(af, na.rm = TRUE),
          ) -> .variant_somatic

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(barcode) |>
          dplyr::summarize(
            somatic_af_cell = mean(af, na.rm = TRUE),
          ) -> .variant_somatic_cell

        .haplo_violin |>
          dplyr::filter(variant %in% .somatic) |>
          dplyr::group_by(cluster) |>
          dplyr::summarize(
            somatic_af_cluster = mean(af, na.rm = TRUE),
          ) -> .cluster_somatic_cluster

        cli::cli_alert_success("Completed processing {i}/{total}")

        tibble::tibble(
          haplo_af = mean(.variant_non_editing$mean_af, na.rm = TRUE),
          somatic_af = mean(.variant_somatic$mean_af, na.rm = TRUE),
          haplo_af_cell = list(.barcode_non_editing_cell),
          somatic_af_cell = list(.variant_somatic_cell),
          haplo_af_cluster = list(.cluster_non_editing_cluster),
          somatic_af_cluster = list(.cluster_somatic_cluster),
        )
      },
      .haplo_violin = haplo_violin,
      .somatic_variant = somatic_variant,
      i = seq_along(haplo_violin),
      total = length(haplo_violin),
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) -> for_hetero_af

for_hetero_af |>
  dplyr::select(gseid, srrid, chemistry, mean_heteroplasmy_count) |>
  tidyr::unnest(cols = mean_heteroplasmy_count) |>
  dplyr::left_join(
    gse_dataset_metadata_full |>
      dplyr::select(srrid, Age, Age_new, Age_group, disease),
    by = c("srrid" = "srrid")
  ) -> for_hetero_af_forplot


# ! somatic --------------------------------------------------------------------

outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz"
celltypes <- c("B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T")
\() {
  for_hetero_af_forplot |>
    dplyr::filter(Age_group != "Unknown") |>
    dplyr::filter(!is.na(somatic_af)) |>
    dplyr::filter(!is.na(Age_new)) |>
    tidyr::unnest(cols = somatic_af_cluster) |>
    dplyr::mutate(
      cluster = factor(cluster, levels = celltypes)
    ) -> for_hetero_af_forplot_cluster

  for_hetero_af_forplot_cluster |>
    dplyr::group_by(Age_group) |>
    dplyr::summarise(
      mean_somatic_af = mean(somatic_af_cluster, na.rm = TRUE)
    ) -> cluster_mean

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
          fill = RColorBrewer::brewer.pal(8, "Set2")
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      ),
      scales = "free_y",
    ) +
    # geom_violin(
    #   aes(
    #     y = somatic_af_cluster,
    #     fill = mean_somatic_af
    #   ),
    #   alpha = 0.5,
    #   size = 1,
    #   color = NA,
    #   show.legend = FALSE
    # ) +
    # scale_fill_gradient2(
    #   name = "AF",
    #   low = "white",
    #   mid = "red",
    #   high = "#3B0049",
    #   midpoint = 0.5,
    # ) +
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
    ) -> p_age_somatic_af_cluster

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
    ) -> for_hetero_af_forplot_cluster

  for_hetero_af_forplot_cluster |>
    dplyr::group_by(Age_group) |>
    dplyr::summarise(
      mean_haplo_af = mean(haplo_af_cluster, na.rm = TRUE)
    ) -> cluster_mean

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
          fill = RColorBrewer::brewer.pal(8, "Set2")
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        by_layer_y = FALSE,
      ),
      scales = "free_y",
    ) +
    # geom_violin(
    #   aes(
    #     y = haplo_af_cluster,
    #     fill = mean_haplo_af
    #   ),
    #   alpha = 0.5,
    #   size = 1,
    #   color = NA,
    #   show.legend = FALSE
    # ) +
    # scale_fill_gradient2(
    #   name = "AF",
    #   low = "white",
    #   mid = "red",
    #   high = "#3B0049",
    #   midpoint = 0.5,
    # ) +
    # ggnewscale::new_scale_fill() +
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
    ) -> p_age_haplot_af_cluster

  ggsave(
    file.path(outdir, "p_age_haplot_af_cluster.pdf"),
    p_age_haplot_af_cluster,
    width = 13,
    height = 7
  )
}


# ! disease and mean heteroplasmy frequency --------------------------------------------------------------------

for_hetero_af_forplot |>
  dplyr::mutate(
    somatic_af = ifelse(
      is.na(somatic_af),
      0,
      somatic_af
    )
  ) |>
  # dplyr::filter(
  #   !is.na(somatic_af)
  # ) |>
  dplyr::filter(
    disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19")
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Alzheimer's Disease",
        "Healthy",
        "COVID-19"
      )
    )
  ) -> for_mhc_disease_boxplot

my_comparisons <- list(
  c("Alzheimer's Disease", "Healthy"),
  c("Alzheimer's Disease", "COVID-19"),
  c("Healthy", "COVID-19")
)

for_mhc_disease_boxplot |>
  ggpubr::ggboxplot(
    x = "disease",
    y = "somatic_af",
    color = "disease",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    add = "jitter"
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),
    legend.position = "none"
  ) +
  labs(
    y = "Mean somatic variant AF",
    title = "All samples"
  ) -> p1
for_mhc_disease_boxplot |>
  dplyr::filter(
    somatic_af > 0
  ) |>
  ggpubr::ggboxplot(
    x = "disease",
    y = "somatic_af",
    color = "disease",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    add = "jitter"
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),
    legend.position = "none"
  ) +
  labs(
    y = "Mean somatic variant AF",
    title = "Samples with # of somatic variants > 0"
  ) -> p2

wrap_plots(list(p1, p2), ncol = 2) -> p_combined
p_combined

ggsave(
  file.path(outdir, "mhc_disease_boxplot_combined.pdf"),
  p_combined,
  width = 10,
  height = 5
)

# ! disease and mean heteroplasmy frequency --------------------------------------------------------------------

gse_dataset_metadata_full |>
  # dplyr::filter(!gseid %in% gseids_tobe_excluded) |>
  # dplyr::filter(gseid != "GSE220189") |>
  dplyr::filter(
    disease %in% c("Alzheimer's Disease", "Healthy", "COVID-19")
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Alzheimer's Disease",
        "Healthy",
        "COVID-19"
      )
    )
  ) -> for_count_disease_boxplot

for_count_disease_boxplot |>
  ggpubr::ggboxplot(
    x = "disease",
    y = "# of somatic variants",
    color = "disease",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    add = "jitter"
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),
    legend.position = "none"
  ) +
  labs(
    y = "# of somatic variants",
    title = "All samples"
  ) -> p1_count

for_count_disease_boxplot |>
  dplyr::filter(
    `# of somatic variants` > 0
  ) |>
  ggpubr::ggboxplot(
    x = "disease",
    y = "# of somatic variants",
    color = "disease",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    add = "jitter"
  ) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),
    legend.position = "none"
  ) +
  labs(
    y = "# of somatic variants",
    title = "Samples with # of somatic variants > 0"
  ) -> p2_count


wrap_plots(list(p1_count, p2_count), ncol = 2) -> p_combined_count
p_combined_count

ggsave(
  file.path(outdir, "somatic_variant_disease_boxplot_combined.pdf"),
  p_combined_count,
  width = 10,
  height = 5
)
# ! somatic --------------------------------------------------------------------

for_hetero_af_forplot |>
  dplyr::mutate(
    somatic_af = ifelse(
      is.na(somatic_af),
      0,
      somatic_af
    )
  ) |>
  dplyr::mutate(
    Age_group = purrr::map_chr(
      .x = Age,
      .f = \(.Age) {
        .x <- as.integer(.Age)
        if (is.na(.x)) {
          return(.x)
        }
        if (.x < 20) {
          return("<20")
        } else if (.x >= 80) {
          return(">80")
        } else {
          # For ages 20-79, use 10-year intervals
          lower <- 10 * floor(.x / 10)
          upper <- lower + 9
          return(paste0(lower, "-", upper))
        }
      }
    )
  ) |>
  dplyr::mutate(
    Age_group = ifelse(is.na(Age_group), "Unknown", Age_group)
  ) |>
  dplyr::mutate(
    Age_group = factor(
      Age_group,
      levels = c(
        "<20",
        "20-29",
        "30-39",
        "40-49",
        "50-59",
        "60-69",
        "70-79",
        ">80",
        "Unknown"
      )
    )
  ) |>
  dplyr::filter(Age_group != "Unknown") |>
  # dplyr::filter(somatic_af > 0) |>
  dplyr::filter(!is.na(Age_new)) -> for_age_plot

for_age_plot |>
  dplyr::group_by(Age_group) |>
  dplyr::summarise(
    mean_somatic_af = mean(somatic_af, na.rm = TRUE)
  ) -> age_mean

for_age_plot |>
  dplyr::left_join(
    age_mean,
    by = c("Age_group" = "Age_group")
  ) |>
  ggplot(aes(x = Age_group)) +
  geom_hline(
    yintercept = 0.3,
    linetype = 21,
    color = "red",
  ) +
  # ggh4x::facet_wrap2(
  #   ~chemistry,
  #   ncol = 1,
  #   strip.position = "right",
  #   strip = ggh4x::strip_themed(
  #     background_y = ggh4x::elem_list_rect(
  #       fill = RColorBrewer::brewer.pal(8, "Set2")
  #     ),
  #     text_y = ggh4x::elem_list_text(
  #       colour = "white",
  #       face = c("bold")
  #     ),
  #     by_layer_y = FALSE,
  #   ),
  #   scales = "free_y",
  # ) +
  # geom_violin(
  #   aes(
  #     y = somatic_af,
  #     fill = mean_somatic_af
  #   ),
  #   alpha = 0.5,
  #   size = 1,
  #   color = NA,
  #   show.legend = FALSE
  # ) +
  # scale_fill_gradient2(
  #   name = "AF",
  #   low = "white",
  #   mid = "red",
  #   high = "#3B0049",
  #   midpoint = 0.5,
  # ) +
  # ggnewscale::new_scale_fill() +
  ggbeeswarm::geom_quasirandom(
    aes(
      x = Age_group,
      y = somatic_af,
      color = somatic_af
    ),
    size = 1,
    dodge.width = .3,
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
    labels = c(0, 0.25, 0.3, 0.5, 0.75, 1),
    breaks = c(0, 0.25, 0.3, 0.5, 0.75, 1),
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
  ) -> p_age_somatic_af
p_age_somatic_af

ggsave(
  file.path(outdir, "p_age_somatic_af.pdf"),
  p_age_somatic_af,
  width = 13,
  height = 7
)


# ! age and somatic and disease--------------------------------------------------------------------

gse_dataset_metadata_full |>
  # dplyr::filter(!gseid %in% gseids_tobe_excluded) |>
  dplyr::filter(!is.na(Age_new)) |>
  dplyr::mutate(
    disease = dplyr::case_when(
      disease %in%
        c("Alzheimer's Disease", "Healthy", "COVID-19", "Unknown") ~ disease,
      TRUE ~ "Other"
    )
  ) |>
  dplyr::mutate(
    disease = factor(
      disease,
      levels = c(
        "Healthy",
        "Alzheimer's Disease",
        "COVID-19",
        "Unknown",
        "Other"
      )
    )
  ) |>
  ggplot(aes(
    x = Age_new,
    y = `# of somatic variants`,
    color = disease
  )) +
  geom_point() +
  ggsci::scale_color_aaas(
    name = "Disease"
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
    # axis.title.x = element_blank(),
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
    ),
    axis.title = element_text(
      size = 16,
      color = "black"
    )
  ) +
  labs(
    x = "Age",
    y = "# of somatic variants"
  ) -> p_age_somatic_af_disease
p_age_somatic_af_disease

ggsave(
  file.path(outdir, "p_age_somatic_af_disease.pdf"),
  p_age_somatic_af_disease,
  width = 12,
  height = 7
)

# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
