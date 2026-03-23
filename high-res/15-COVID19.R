#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 01:16:11
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Source ---------------------------------------------------------------------
source("high-res/00-colors.R")
source("high-res/plot_mtdna.R")
# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})


outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))

ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
variant_annotation <- import(outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx")

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------
fn_plot_disease_chemistry <- function(meta) {
  meta |>
    count(disease, Chemistry) |>
    ggplot(aes(x = disease, y = n, fill = Chemistry)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_chemistry) +
    theme_bw() +
    labs(x = "Disease", y = "Count", fill = "Chemistry") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    )
}
fn_plot_ad_variant_bar <- function(df_wide, annotation = NULL) {
  color_disease_category <- c(
    "AD / PD" = "#ff1e05",
    "LHON" = "#b300ff",
    "Metabolic / Diabetes" = "#fe7702",
    "Exercise / Altitude" = "#00d0ff",
    "Neurological" = "#0099ff",
    "Other" = "grey55"
  )

  # df_wide: variant_nsamples_wide (has Healthy, AD, total, start columns)
  # long format for geom_col
  df_wide |>
    dplyr::arrange(-total) |>
    dplyr::mutate(
      variant = factor(variant, levels = variant)
    ) -> df_sorted

  df_sorted |>
    tidyr::pivot_longer(
      cols = c(Healthy, `COVID-19`),
      names_to = "disease",
      values_to = "nsamples"
    ) |>
    dplyr::mutate(
      disease = factor(disease, levels = c("Healthy", "COVID-19"))
    ) -> df_long

  fn_xy_breaks_limits(df_sorted$total, step = 10) -> xyb

  p <- df_long |>
    ggplot(aes(x = variant, y = nsamples)) +
    geom_col(aes(fill = disease)) +
    scale_fill_manual(name = "Disease", values = color_disease) +
    geom_hline(
      yintercept = c(2, xyb$breaks),
      color = "grey80",
      linetype = "dashed"
    )

  if (!is.null(annotation)) {
    annotation |>
      filter(
        prediction_class %in%
          c("Likely pathogenic", "Pathogenic") |
          !is.na(Disease)
      ) -> annotation
    p <- p +
      ggrepel::geom_label_repel(
        data = df_sorted |>
          dplyr::inner_join(
            annotation |> dplyr::select(variant, label, disease_category),
            by = "variant"
          ),
        aes(
          x = variant,
          y = total,
          label = label,
          color = disease_category,
          segment.colour = after_scale(colour)
        ),
        inherit.aes = FALSE,
        size = 2.5,
        lineheight = 0.9,
        label.size = 0.2,
        label.padding = unit(0.2, "lines"),
        label.r = unit(0.1, "lines"),
        box.padding = 0.4,
        point.padding = 0.3,
        min.segment.length = 0,
        segment.size = 0.3,
        segment.curvature = -0.2,
        direction = "both",
        nudge_y = max(xyb$breaks) * 0.15,
        max.overlaps = Inf,
        fill = alpha("white", 0.85)
      ) +
      scale_color_manual(
        name = "Variant annotation",
        values = color_disease_category
      )
  }

  p +
    scale_x_discrete(
      expand = expansion(add = c(3, 1)),
      name = "Variant (sorted by total samples)"
    ) +
    scale_y_continuous(
      breaks = xyb$breaks,
      limits = xyb$limits,
      labels = scales::comma_format(),
      expand = expansion(mult = c(0.01, 0.05)),
      name = "# of Samples"
    ) +
    guides(
      fill = guide_legend(ncol = 1, order = 1),
      color = guide_legend(ncol = 1, order = 2)
    ) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = c(0.7, 0.8),
      legend.box = "horizontal",
      legend.key = element_blank(),
      axis.title.y = element_text(size = 16, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      legend.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 16, colour = "black"),
      strip.background = element_blank(),
      plot.title = element_text(
        size = 16,
        color = "black",
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(title = "AD Variant Distribution (ranked by frequency)")
}
# Main --------------------------------------------------------------------

METAFULL |>
  filter(
    disease %in% c("Healthy", "COVID-19"),
    # Chemistry == "SC5P-PE"
  ) |>
  select(gseid, srrid, Chemistry, disease) -> disease_meta

disease_meta |>
  fn_plot_disease_chemistry() +
  labs(
    title = "Distribution of Disease and Chemistry"
  ) -> plot_covid19_disease_chemistry

saveplot(
  filename = outdirnotuse /
    "COVID19" /
    "COVID19-DISEASE-CHEMISTRY-DISTRIBUTION.pdf",
  plot = plot_covid19_disease_chemistry,
  width = 6,
  height = 6
)


disease_meta |>
  dplyr::left_join(ALLVARIANTS, by = c("gseid", "srrid")) |>
  dplyr::select(-c(Chemistry, Haplogroup, Verbose_haplogroup)) |>
  dplyr::mutate(
    disease = factor(disease, levels = c("Healthy", "COVID-19"))
  ) -> meta_af


meta_af |>
  tidyr::nest(.by = c(variant, disease)) |>
  dplyr::mutate(nsamples = purrr::map_int(data, nrow)) |>
  dplyr::select(-data) |>
  tidyr::pivot_wider(
    names_from = disease,
    values_from = nsamples,
    values_fill = 0
  ) |>
  dplyr::mutate(
    total = `Healthy` + `COVID-19`,
    start = as.integer(gsub("(\\d+)[A-Z]>[A-Z]", "\\1", variant))
  ) |>
  dplyr::arrange(-total) -> variant_nsamples_wide


# -- Variant annotation labels --
load_pkg(stringr)
variant_annotation |>
  dplyr::filter(variant %in% meta_af$variant) |>
  # filter(
  #   prediction_class %in% c("pathogenic", "likely_pathogenic") | !is.na(Disease)
  # ) |>
  dplyr::mutate(
    label = glue::glue("{variant}\n{Disease}\n{prediction_class}"),
    disease_category = dplyr::case_when(
      str_detect(
        Disease,
        regex("\\bAD\\b|\\bPD\\b|Alzheimer|Parkinson", ignore_case = TRUE)
      ) ~ "AD / PD",
      str_detect(Disease, regex("LHON", ignore_case = TRUE)) ~ "LHON",
      str_detect(
        Disease,
        regex("T2D|diabetes|metabolic", ignore_case = TRUE)
      ) ~ "Metabolic / Diabetes",
      str_detect(
        Disease,
        regex("altitude|VO2|exercise|EXIT|cyclic vomiting", ignore_case = TRUE)
      ) ~ "Exercise / Altitude",
      str_detect(
        Disease,
        regex(
          "DEAF|SNHL|hearing|dystonia|encephalomyopathy|Mitochondrial|neuropathy|Respiratory Chain",
          ignore_case = TRUE
        )
      ) ~ "Neurological",
      TRUE ~ "Other"
    )
  ) -> disease_variant_annotation

# -- Export intermediate data --
variant_nsamples_wide |>
  export(outdirnotuse / "COVID19" / "COVID19-variant-sample-counts.xlsx")
disease_meta |> export(outdirnotuse / "COVID19" / "COVID19-variant-af.xlsx")
disease_meta |> export(outdirnotuse / "COVID19" / "COVID19-variant-af.qs")


# Plot 2: sorted bar chart (previous style)
fn_plot_ad_variant_bar(
  variant_nsamples_wide,
  disease_variant_annotation
) -> p_disease_variant_bar


saveplot(
  plot = p_disease_variant_bar,
  filename = outdirnotuse /
    "COVID19" /
    "COVID19-variant-distribution-ranked.pdf",
  width = 17,
  height = 8
)


variant_nsamples_wide |>
  dplyr::arrange(-total) |>
  dplyr::mutate(
    variant = factor(variant, levels = variant)
  ) |>
  dplyr::left_join(
    disease_variant_annotation,
    by = "variant"
  ) -> variant_nsamples_wide_annotated
variant_nsamples_wide_annotated |>
  export(
    outdirnotuse / "COVID19" / "COVID19-variant-sample-counts-annotated.xlsx"
  )

variant_nsamples_wide_annotated |>
  filter(
    `COVID-19` > 5,
    `Healthy` > 5
  ) |>
  count(Disease, prediction_class)
variant_nsamples_wide_annotated |>
  filter(
    `COVID-19` > 5,
    `Healthy` > 5
  ) |>
  filter(!is.na(Disease)) |>
  glimpse()


variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  )


variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  )


variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  ) |>
  filter(
    `COVID-19` > 1
  )

variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  ) |>
  filter(
    `COVID-19` > 1
  ) |>
  filter(!is.na(Disease)) |>
  glimpse()

variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  ) |>
  filter(
    `COVID-19` > 1
  ) |>
  count(Disease, prediction_class)
#

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
