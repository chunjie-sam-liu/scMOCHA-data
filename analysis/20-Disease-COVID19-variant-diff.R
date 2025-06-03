#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-09 16:26:12
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
outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/disease"

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

gse_dataset_metadata_full <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)

gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)

all_variant <- import("/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.qs") |>
  dplyr::select(variant, issomatic)

all_hetero_af_cluster <- import(
  file.path(cleandatadir, "all_hetero_af.cluster.csv"),
) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, barcode, num_variants),
    names_to = "variant",
    values_to = "af"
  )

all_hetero_af_bulk <- import(
  file.path(cleandatadir, "all_hetero_af.bulk.csv"),
) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, barcode, num_variants),
    names_to = "variant",
    values_to = "af"
  )

all_hetero_af_cluster |>
  dplyr::bind_rows(
    all_hetero_af_bulk
  ) ->
all_hetero_af_cluster_bulk

# body --------------------------------------------------------------------


gse_dataset_metadata_full |>
  dplyr::filter(
    disease %in% c("Healthy", "COVID-19")
  ) |>
  dplyr::select(
    gseid, srrid, Chemistry, disease
  ) ->
admeta

admeta |>
  dplyr::count(
    disease, Chemistry
  )

admeta |>
  dplyr::filter(
    Chemistry != "SC5P-PE"
  ) ->
admeta_sc5p


admeta_sc5p |>
  dplyr::left_join(
    gse_data,
    by = c("gseid", "srrid")
  ) ->
admeta_sc5p_variant

admeta_sc5p_variant |>
  dplyr::mutate(
    variant_type = purrr::map(
      .x = anno,
      .f = function(.x) {
        .x |>
          dplyr::mutate(
            variant = glue::glue("{Position}{Ref}>{Alt}")
          ) |>
          dplyr::select(variant, ntchange) |>
          dplyr::left_join(
            all_variant,
            by = "variant"
          )
      }
    )
  ) |>
  dplyr::select(
    gseid, srrid, Chemistry, disease, variant_type
  ) ->
admeta_sc5p_variant_type




# ! compare variant between AD and healthy --------------------------------------------------------------------

admeta_sc5p_variant_type |>
  tidyr::unnest(cols = variant_type) |>
  dplyr::filter(
    issomatic == "heteroplasmic"
  ) |>
  dplyr::select(
    srrid, disease, variant
  ) |>
  dplyr::left_join(
    all_hetero_af_cluster_bulk,
    by = c("srrid", "variant")
  ) ->
admeta_sc5p_variant_type_af


admeta_sc5p_variant_type_af |>
  dplyr::group_by(
    variant, barcode
  ) |>
  tidyr::nest() |>
  dplyr::ungroup() |>
  dplyr::mutate(
    t = purrr::map(
      .x = data,
      .f = \(.x) {
        # .x <- a$data[[1]]
        .x |>
          dplyr::count(disease) |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n
          ) ->
        .xx
        tryCatch(
          expr = {
            t.test(
              af ~ disease,
              data = .x,
              var.equal = TRUE
            ) |>
              broom::tidy() |>
              dplyr::select(
                estimate, estimate1, estimate2, p.value, conf.low, conf.high
              ) |>
              dplyr::bind_cols(
                .xx
              )
          },
          error = function(e) {
            message("Error: ", conditionMessage(e))
            return(NULL)
          }
        )
      }
    )
  ) ->
admeta_sc5p_variant_type_af_ttest


admeta_sc5p_variant_type_af_ttest |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = t) |>
  dplyr::filter(p.value < 0.2) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    est = abs(estimate),
  ) |>
  dplyr::mutate(
    rank = plog10p * est,
  ) |>
  dplyr::arrange(
    desc(rank)
  ) |>
  dplyr::rename(
    ad = "COVID-19",
  ) |>
  dplyr::filter(
    ad >= 10,
    Healthy >= 10
  ) ->
admeta_sc5p_variant_type_af_ttest_rank

admeta_sc5p_variant_type_af_ttest_rank |>
  dplyr::filter(estimate > 0.03) |>
  dplyr::arrange(variant) ->
top5_variant


topvariants <- c("1397T>A", "1670A>G", "3173G>A", "3176A>T", "3178T>A")


source("/home/liuc9/github/scMOCHA-data/analysis/00-colors.R")

library(ggh4x)
library(ggbeeswarm)
library(ggnewscale)

color_celltype_bulk <- c(
  "Pseudo-bulk" = "red",
  color_celltype
)

admeta_sc5p_variant_type_af_ttest |>
  dplyr::select(-data) |>
  tidyr::unnest(cols = t) |>
  dplyr::rename(
    ad = "COVID-19",
  ) |>
  # dplyr::filter(
  #   ad >= 5,
  #   Healthy >= 5
  # ) |>
  dplyr::mutate(
    log10p = -log10(p.value),
  ) |>
  dplyr::mutate(
    barcode = gsub(barcode, pattern = "_", replacement = " "),
    barcode = ifelse(barcode == "bulk", "Pseudo-bulk", barcode),
  ) |>
  dplyr::mutate(
    barcode = factor(barcode, levels = names(color_celltype_bulk)),
  ) |>
  dplyr::mutate(
    point_color = ifelse(
      p.value < 0.05 & abs(estimate) > 0.03,
      "red",
      "black"
    ),
  ) |>
  dplyr::mutate(
    label = ifelse(
      p.value < 0.05 & abs(estimate) > 0.03,
      variant,
      NA
    ),
  ) ->
forplot_test

# forplot_test |>
#   dplyr::filter(variant == "1670A>G")

forplot_test |>
  ggplot(aes(
    x = estimate,
    y = log10p,
  )) +
  geom_point(
    aes(
      color = point_color,
    ),
  ) +
  scale_color_identity() +
  ggrepel::geom_text_repel(
    aes(
      label = label,
    ),
    size = 3,
    show.legend = FALSE,
    # nudge_x = 0.1,
    # nudge_y = 0.1,
    # segment.size = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 20,
    color = "red"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 20,
    color = "red"
  ) +
  ggh4x::facet_grid2(
    ~barcode,
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = color_celltype_bulk,
        color = NA
      ),
      text_x = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      )
    )
  ) +
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    # panel.background = element_blank(),
    panel.background = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    # axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
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
    y = "-log10(p-value)",
    x = "Effect Size (COVID-19 - Healthy)",
  ) ->
p_variant_boxplot_af_sc_ttest


\(){
  admeta_sc5p_variant_type_af |>
    dplyr::filter(variant %in% topvariants) |>
    dplyr::mutate(
      barcode = gsub(barcode, pattern = "_", replacement = " "),
      barcode = ifelse(barcode == "bulk", "Pseudo-bulk", barcode),
    ) |>
    dplyr::mutate(
      barcode = factor(barcode, levels = names(color_celltype_bulk)),
    ) ->
  forplot_

  forplot_ |>
    dplyr::group_by(disease, barcode, variant) |>
    dplyr::summarise(
      mean_cluster_variant_af = mean(af, na.rm = T),
    ) ->
  disease_barcode_variatn_af

  forplot_ |>
    dplyr::left_join(
      disease_barcode_variatn_af,
      by = c("disease", "barcode", "variant")
    ) |>
    ggplot(aes(x = disease)) +
    ggh4x::facet_grid2(
      variant ~ barcode,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = color_celltype_bulk,
          color = NA
        ),
        text_x = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        ),
        background_y = ggh4x::elem_list_rect(
          fill = "black",
          color = NA
        ),
        text_y = ggh4x::elem_list_text(
          colour = "white",
          face = c("bold")
        )
      )
    ) +
    geom_violin(
      aes(
        y = af,
        fill = disease
      ),
      alpha = 0.7,
      size = 1,
      color = NA
    ) +
    scale_fill_manual(
      values = color_disease[c("Alzheimer's Disease", "Healthy")],
      name = "Disease",
      labels = c("AD", "Healthy")
    ) +
    ggbeeswarm::geom_quasirandom(
      aes(
        y = af,
        color = disease
      ),
      size = 1,
      dodge.width = .75,
      alpha = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = color_disease[c("Alzheimer's Disease", "Healthy")],
      name = "Disease",
      labels = c("AD", "Healthy")
    ) +
    ggsignif::geom_signif(
      aes(
        y = af,
      ),
      comparisons = list(
        c("Alzheimer's Disease", "Healthy")
      ),
      y_position = 0.8
    ) +
    theme(
      plot.margin = margin(t = 0, b = 0, unit = "cm"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      # axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "top",
      legend.key = element_blank(),
      axis.title.y = element_text(color = "black", size = 16),
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
      y = "Allele Frequency",
    ) ->
  p_variant_boxplot_af

  ggsave(
    filename = file.path(outdir, "covid19-variant_boxplot.pdf"),
    plot = p_variant_boxplot_af,
    width = 15,
    height = 8,
    device = cairo_pdf
  )
}


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
