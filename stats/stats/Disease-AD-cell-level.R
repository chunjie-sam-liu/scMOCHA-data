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
outdir <- "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/heteroplasmic"

cleandatadir <- "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"

gse_dataset_metadata_full <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_dataset_metadata_full.rds"
)

gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.rds"
)

all_variant <- readr::read_rds("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/all_variant.rds") |>
  dplyr::select(variant, issomatic)

all_heteroplasmic_af <- data.table::fread(
  file.path(cleandatadir, "all_hetero_af.cluster.csv"),
  header = TRUE,
  sep = ",",
) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, barcode, num_variants),
    names_to = "variant",
    values_to = "af"
  )

# body --------------------------------------------------------------------


gse_dataset_metadata_full |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease")
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
    # Chemistry %in% c("SC5P-PE", "SC5P-R2")
    Chemistry == "SC5P-PE"
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



# ! variant venn --------------------------------------------------------------------

\(){
  admeta_sc5p_variant_type |>
    tidyr::unnest(cols = variant_type) |>
    dplyr::filter(
      issomatic == "heteroplasmic"
    ) |>
    dplyr::select(
      srrid, disease, variant
    ) |>
    dplyr::group_by(
      variant
    ) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      m = purrr::map(
        .x = data,
        .f = function(.x) {
          .x |>
            dplyr::group_by(disease) |>
            dplyr::count() |>
            dplyr::ungroup() |>
            tidyr::pivot_wider(
              names_from = disease,
              values_from = n
            )
        }
      )
    ) |>
    tidyr::unnest(cols = m) ->
  admeta_sc5p_variant_type_count

  admeta_sc5p_variant_type_count |>
    dplyr::select(-data) |>
    dplyr::arrange(
      Healthy
    ) |>
    dplyr::filter(
      variant == "3173G>A"
    )
}


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
    all_heteroplasmic_af,
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
  dplyr::filter(p.value < 0.05) |>
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
  dplyr::filter(
    `Alzheimer's Disease` > 10,
    Healthy > 10
  ) |>
  dplyr::filter(barcode == "Mono") ->
admeta_sc5p_variant_type_af_ttest_rank

admeta_sc5p_variant_type_af_ttest_rank |>
  dplyr::filter(estimate > 0.03) |>
  dplyr::arrange(variant) ->
top5_variant

topvariants <- c("1397T>A", "1670A>G", "3173G>A", "3176A>T", "3178T>A")


source("/home/liuc9/github/scMOCHA-data/stats/stats/00-colors.R")
color_chemistry

admeta_sc5p_variant_type_af |>
  dplyr::filter(variant %in% topvariants) |>
  ggplot(aes(
    x = disease, y = af, fill = disease
  )) +
  geom_boxplot(outlier.color = NA, width = 0.5) +
  geom_jitter(
    size = 0.5,
    alpha = 0.5,
    width = 0.2
  ) +
  scale_fill_manual(
    values = color_disease[c("Alzheimer's Disease", "Healthy")],
    name = "Disease",
    labels = c("AD", "Healthy")
  ) +
  ggsignif::geom_signif(
    comparisons = list(
      c("Alzheimer's Disease", "Healthy")
    ),
    y_position = 0.8
  ) +
  ggh4x::facet_grid2(
    variant ~ barcode,
    # switch = "x",
    strip = ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = color_celltype
      ),
      text_x = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      ),
      background_y = ggh4x::elem_list_rect(
        fill = "#AAAAAA"
      ),
      text_y = ggh4x::elem_list_text(
        colour = "white",
        face = c("bold")
      )
    )
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
    legend.position = "right",
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
  filename = file.path(outdir, "variant_boxplot_af.pdf"),
  plot = p_variant_boxplot_af,
  width = 15,
  height = 8,
  device = cairo_pdf
)




admeta_sc5p_variant_type_af |>
  dplyr::filter(
    variant == admeta_sc5p_variant_type_af_ttest_rank$variant[3],
    barcode == admeta_sc5p_variant_type_af_ttest_rank$barcode[3]
  ) |>
  ggstatsplot::ggbetweenstats(
    data = _,
    x = disease,
    y = af,
    pairwise.display = "p-value",
    pairwise.comparisons = TRUE,
    p.adjust.method = "fdr",
    p.adjust.display = TRUE,
    p.value.label = "p.adj",
    p.value.label.size = 3.5,
    p.value.label.color = "black",
    p.value.label.position = c(0.5, 0.95),
    p.value.label.nudge_y = 0.05,
    ggplot.component = list(
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    )
  )

pcc <- readr::read_tsv(file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv") |>
  dplyr::arrange(cancer_types)
admeta_sc5p_variant_type_af |>
  dplyr::filter(
    variant == "1670A>G"
  ) |>
  dplyr::filter(
    variant %in% top5_variant$variant
  ) |>
  ggpubr::ggboxplot(
    data = _,
    x = "disease",
    y = "af",
    add = "jitter",
    color = "disease",
    palette = "jco",
    facet.by = "barcode",
    nrow = 1,
    ggplot.component = list(
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    )
  ) +
  ggpubr::stat_compare_means()


# footer ------------------------------------------------------------------

# future: :plan(future: :sequential)

# save image --------------------------------------------------------------
