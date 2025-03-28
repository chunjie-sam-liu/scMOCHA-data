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


# for_hetero_af_forplot |>
#   dplyr::filter(Age_group != "Unknown") |>
#   dplyr::filter(!is.na(haplo_af)) |>
#   ggplot(aes(
#     x = Age_group,
#     y = haplo_af,
#   )) +
#   geom_boxplot() +
#   facet_wrap(
#     ~chemistry,
#     ncol = 3
#   )


# for_hetero_af_forplot |>
#   dplyr::filter(Age_group != "Unknown") |>
#   dplyr::filter(!is.na(somatic_af)) |>
#   dplyr::filter(!is.na(Age_new)) |>
#   # dplyr::filter(disease == "Healthy") |>
#   ggplot(aes(
#     x = Age_group,
#     y = haplo_af,
#   )) +
#   geom_boxplot() +
#   facet_wrap(
#     ~chemistry,
#     ncol = 3
#   )

# for_hetero_af_forplot |>
#   dplyr::filter(Age_group != "Unknown") |>
#   dplyr::mutate(
#     somatic_af = ifelse(
#       is.na(somatic_af),
#       0,
#       somatic_af
#     )
#   ) |>
#   dplyr::filter(!is.na(somatic_af)) |>
#   dplyr::filter(!is.na(Age_new)) |>
#   # dplyr::filter(disease == "Healthy") |>
#   ggplot(aes(
#     x = Age_group,
#     y = somatic_af,
#   ))




# ! somatic --------------------------------------------------------------------

outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz"
celltypes <- c("B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T")


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




# ! disease and mha --------------------------------------------------------------------

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
  ) ->
for_mhc_disease_boxplot

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
    y = "Mean somatic variant AF"
  ) ->
p1
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
    y = "Mean somatic variant AF"
  ) ->
p2

wrap_plots(list(p1, p2), ncol = 2) -> p_combined

ggsave(
  file.path(outdir, "mhc_disease_boxplot_combined.pdf"),
  p_combined,
  width = 10,
  height = 5
)

gse_dataset_metadata_full |>
  dplyr::filter(gseid != "GSE220189") |>
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
  ) ->
for_count_disease_boxplot

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
    y = "# of somatic variants"
  ) ->
p1_count

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
    y = "# of somatic variants"
  ) ->
p2_count


wrap_plots(list(p1_count, p2_count), ncol = 2) -> p_combined_count

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
    Age_group = factor(Age_group, levels = c(
      "<20", "20-29", "30-39", "40-49", "50-59",
      "60-69", "70-79", ">80", "Unknown"
    ))
  ) |>
  dplyr::filter(Age_group != "Unknown") |>
  # dplyr::filter(somatic_af > 0) |>
  dplyr::filter(!is.na(Age_new)) ->
for_age_plot

for_age_plot |>
  dplyr::group_by(Age_group) |>
  dplyr::summarise(
    mean_somatic_af = mean(somatic_af, na.rm = TRUE)
  ) ->
age_mean

for_age_plot |>
  dplyr::left_join(
    age_mean,
    by = c("Age_group" = "Age_group")
  ) |>
  ggplot(aes(x = Age_group)) +
  geom_hline(
    yintercept = 0.3,
    linetype = 21,
    color = "black",
  ) +
  # ggh4x::facet_wrap2(
  #   ~chemistry,
  #   ncol = 1,
  #   strip.position = "right",
  #   strip = ggh4x::strip_themed(
  #     background_y = ggh4x::elem_list_rect(
  #       fill = pcc$color
  #     ),
  #     text_y = ggh4x::elem_list_text(
  #       colour = "white",
  #       face = c("bold")
  #     ),
  #     by_layer_y = FALSE,
  #   ),
  #   scales = "free_y",
  # ) +
  geom_violin(
    aes(
      y = somatic_af,
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
p_age_somatic_af

ggsave(
  file.path(outdir, "p_age_somatic_af.pdf"),
  p_age_somatic_af,
  width = 13,
  height = 7
)


# meta --------------------------------------------------------------------

gse_dataset_metadata_full |>
  dplyr::filter(gseid != "GSE220189") |>
  dplyr::select(gseid, srrid, Race, Ethnicity, Gender, Age_group, disease, Chemistry) ->
gse_dataset_metadata_full_selected

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Gender) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Gender_str = glue::glue("{Gender}\n(n={n})")) %>%
  dplyr::mutate(Gender_str = factor(Gender_str, levels = Gender_str)) ->
Gender_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Race) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Race_str = glue::glue("{Race}\n(n={n})")) %>%
  dplyr::mutate(Race_str = factor(Race_str, levels = Race_str)) ->
Race_str


gse_dataset_metadata_full_selected |>
  dplyr::group_by(Ethnicity) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Ethnicity_str = glue::glue("{Ethnicity}\n(n={n})")) %>%
  dplyr::mutate(Ethnicity_str = factor(Ethnicity_str, levels = Ethnicity_str)) ->
Ethnicity_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(disease) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Disease_str = glue::glue("{disease}\n(n={n})")) %>%
  dplyr::mutate(Disease_str = factor(Disease_str, levels = Disease_str)) ->
Disease_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Chemistry) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Chemistry_str = glue::glue("{Chemistry}\n(n={n})")) %>%
  dplyr::mutate(Chemistry_str = factor(Chemistry_str, levels = Chemistry_str)) ->
Chemistry_str

gse_dataset_metadata_full_selected |>
  dplyr::mutate(
    Age = ifelse(Age_group == "Unknown", "Unknown", "Known")
  ) |>
  dplyr::group_by(Age) |>
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Age_str = glue::glue("{Age}\n(n={n})")) %>%
  dplyr::mutate(Age_str = factor(Age_str, levels = Age_str)) ->
Age_str


gse_dataset_metadata_full_selected |>
  dplyr::mutate(
    Age = ifelse(Age_group == "Unknown", "Unknown", "Known")
  ) |>
  dplyr::select(-gseid, -srrid, -Age_group) |>
  dplyr::group_by(Chemistry, Age, Gender, Race, Ethnicity, disease) |>
  dplyr::count() %>%
  dplyr::ungroup() |>
  dplyr::left_join(
    Chemistry_str,
    by = c("Chemistry")
  ) |>
  dplyr::left_join(
    Age_str,
    by = c("Age")
  ) |>
  dplyr::left_join(
    Gender_str,
    by = "Gender"
  ) |>
  dplyr::left_join(
    Race_str,
    by = "Race"
  ) |>
  dplyr::left_join(
    Ethnicity_str,
    by = "Ethnicity"
  ) |>
  dplyr::left_join(
    Disease_str,
    by = "disease"
  ) ->
for_sankey_plot
library(ggalluvial)

chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
ggsci::pal_aaas()(3) |> prismatic::color()
chem_colors <- viridis::viridis_pal(option = "D")(4) |>
  prismatic::color()

for_sankey_plot |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = chem_levels
    )
  ) |>
  ggplot(aes(axis1 = Chemistry_str, axis2 = Age_str, axis3 = Gender_str, axis4 = Race_str, axis5 = Ethnicity_str, axis6 = Disease_str, y = n.x)) +
  ggalluvial::geom_alluvium(aes(fill = Chemistry), width = 1 / 12) +
  ggalluvial::geom_stratum(fill = "white") +
  scale_fill_manual(values = chem_colors) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(
    limits = c("Chemistry_str", "Age_str", "Gender_str", "Race_str", "Ethnicity_str", "Disease_str"),
    labels = gsub("_str", "", c("Chemistry_str", "Age_str", "Gender_str", "Race_str", "Ethnicity_str", "Disease_str")),
    expand = c(0.2, 0.05)
  ) +
  scale_y_continuous(
    expand = expand_scale(mult = c(0.01, 0.02))
  ) +
  # ggsci::scale_fill_aaas() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(title = "Sequencing", nrow = 1)) ->
meta_plot_sankey

ggsave(
  file.path(outdir, "meta_plot_sankey.pdf"),
  meta_plot_sankey,
  width = 15,
  height = 8
)


# ! age and somatic --------------------------------------------------------------------

gse_dataset_metadata_full |>
  dplyr::filter(!is.na(Age_new)) |>
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
  ) ->
p_age_somatic_af_disease

ggsave(
  file.path(outdir, "p_age_somatic_af_disease.pdf"),
  p_age_somatic_af_disease,
  width = 12,
  height = 7
)


gse_dataset_metadata_full |>
  dplyr::filter(!is.na(Age_new)) |>
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
    color = Chemistry
  )) +
  geom_point() +
  scale_color_brewer(palette = "Set1")



# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
