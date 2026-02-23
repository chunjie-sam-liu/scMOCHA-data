#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-09 19:46:38
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


# Source ------------------------------------------------------------------
source("high-res/00-colors.R")
source("high-res/plot_mtdna.R")

# Load data ---------------------------------------------------------------
load_pkg(stringr)
dotenv(".env")
suppressMessages(conflicted::conflict_prefer("filter", "dplyr"))

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))

ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
variant_annotation <- import(outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx")

# Colors ------------------------------------------------------------------
color_disease_category <- c(
  "AD / PD" = "#ff1e05",
  "LHON" = "#b300ff",
  "Metabolic / Diabetes" = "#fe7702",
  "Exercise / Altitude" = "#00d0ff",
  "Neurological" = "#0099ff",
  "Other" = "grey55"
)

# Function ----------------------------------------------------------------
fn_plot_ad_variant_hotspots <- function(df, annotation = NULL) {
  fn_xy_breaks_limits(df$total, step = 10) -> xyb

  df |>
    dplyr::mutate(
      y_healthy_end = Healthy,
      y_ad_start = Healthy,
      y_ad_end = total
    ) -> segs

  p <- segs |>
    ggplot() +
    # Healthy segment: 0 -> Healthy count
    geom_segment(
      aes(x = start, xend = start, y = 0, yend = y_healthy_end),
      color = color_disease[["Healthy"]],
      linewidth = 0.7
    ) +
    # AD segment: Healthy -> total
    geom_segment(
      aes(x = start, xend = start, y = y_ad_start, yend = y_ad_end),
      color = color_disease[["Alzheimer's Disease"]],
      linewidth = 0.7
    ) +
    # Top point sized by total
    geom_point(
      aes(x = start, y = total, size = total),
      color = "grey30",
      shape = 21,
      fill = "grey30",
      stroke = 0.5
    ) +
    # Dummy layer to drive disease fill legend
    geom_point(
      data = data.frame(
        disease = factor(
          c("Healthy", "Alzheimer's Disease"),
          levels = c("Healthy", "Alzheimer's Disease")
        )
      ),
      aes(x = -Inf, y = -Inf, fill = disease),
      shape = 22,
      size = 4,
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = c(2, xyb$breaks),
      color = "grey80",
      linetype = "dashed"
    )

  if (!is.null(annotation)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = df |>
          dplyr::inner_join(
            annotation |> dplyr::select(variant, label, disease_category),
            by = "variant"
          ),
        aes(
          x = start,
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
      )
  }

  p +
    scale_fill_manual(name = "Disease", values = color_disease) +
    scale_color_manual(
      name = "Variant annotation",
      values = color_disease_category
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.01)),
      limits = c(1, 16569),
      breaks = c(seq(0, 17000, 1000), 16569),
      labels = c(seq(0, 17000, 1000), 16569)
    ) +
    scale_y_continuous(
      breaks = xyb$breaks,
      limits = xyb$limits,
      labels = scales::comma_format(),
      expand = expansion(mult = c(0.01, 0.05)),
      name = "# of Samples"
    ) +
    guides(
      size = "none",
      fill = guide_legend(ncol = 1, order = 1),
      color = guide_legend(ncol = 1, order = 2)
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
    labs(title = "mtDNA Variant Hotspots - Healthy vs. Alzheimer's Disease")
}

fn_plot_ad_variant_bar <- function(df_wide, annotation = NULL) {
  # df_wide: variant_nsamples_wide (has Healthy, AD, total, start columns)
  # long format for geom_col
  df_wide |>
    dplyr::arrange(-total) |>
    dplyr::mutate(
      variant = factor(variant, levels = variant)
    ) -> df_sorted

  df_sorted |>
    tidyr::pivot_longer(
      cols = c(Healthy, `Alzheimer's Disease`),
      names_to = "disease",
      values_to = "nsamples"
    ) |>
    dplyr::mutate(
      disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
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

# Body --------------------------------------------------------------------

# -- AD metadata --
METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease"),
    Chemistry == "SC5P-PE"
  ) |>
  dplyr::select(gseid, srrid, Chemistry, disease) -> admeta

admeta |>
  dplyr::left_join(ALLVARIANTS, by = c("gseid", "srrid")) |>
  dplyr::select(-c(Chemistry, Haplogroup, Verbose_haplogroup)) |>
  dplyr::mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta_af

# -- Variant counts per disease with mtDNA start position --
admeta_af |>
  tidyr::nest(.by = c(variant, disease)) |>
  dplyr::mutate(nsamples = purrr::map_int(data, nrow)) |>
  dplyr::select(-data) |>
  tidyr::pivot_wider(
    names_from = disease,
    values_from = nsamples,
    values_fill = 0
  ) |>
  dplyr::mutate(
    total = `Healthy` + `Alzheimer's Disease`,
    start = as.integer(gsub("(\\d+)[A-Z]>[A-Z]", "\\1", variant))
  ) |>
  dplyr::arrange(-total) -> variant_nsamples_wide

# -- Variant annotation labels --
variant_annotation |>
  dplyr::filter(variant %in% admeta_af$variant) |>
  filter(
    prediction_class %in% c("pathogenic", "likely_pathogenic") | !is.na(Disease)
  ) |>
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
  ) -> ad_variant_annotation

# -- Export intermediate data --
variant_nsamples_wide |>
  export(outdirnotuse / "AD" / "AD-variant-sample-counts.xlsx")
admeta_af |> export(outdirnotuse / "AD" / "AD-variant-af.xlsx")
admeta_af |> export(outdirnotuse / "AD" / "AD-variant-af.qs")

# -- Plot & save --
# Plot 1: mtDNA coordinate lollipop
fn_plot_ad_variant_hotspots(
  variant_nsamples_wide,
  ad_variant_annotation
) -> p_ad_variant_dis

wrap_plots(
  p_ad_variant_dis,
  plot_spacer(),
  fn_plot_mtdna(),
  ncol = 1,
  heights = c(15, -0.7, 1)
) -> p_ad_variant_dis_wrapped

ggsave(
  p_ad_variant_dis_wrapped,
  filename = outdirnotuse / "AD" / "AD-variant-distribution-mtdna-coord.pdf",
  width = 17,
  height = 8
)

# Plot 2: sorted bar chart (previous style)
fn_plot_ad_variant_bar(
  variant_nsamples_wide,
  ad_variant_annotation
) -> p_ad_variant_bar

ggsave(
  p_ad_variant_bar,
  filename = outdirnotuse / "AD" / "AD-variant-distribution-ranked.pdf",
  width = 17,
  height = 8
)


# ! don't run below

# T-test: per-variant per-celltype AF difference (Healthy vs AD) ----------
fn_ttest_af <- function(.x) {
  .x |>
    dplyr::count(disease) |>
    tidyr::pivot_wider(names_from = disease, values_from = n) -> counts
  tryCatch(
    t.test(af ~ disease, data = .x, var.equal = TRUE) |>
      broom::tidy() |>
      dplyr::select(
        estimate,
        estimate1,
        estimate2,
        p.value,
        conf.low,
        conf.high
      ) |>
      dplyr::bind_cols(counts),
    error = function(e) NULL
  )
}

admeta_af |>
  tidyr::pivot_longer(
    cols = c(B:Bulk),
    names_to = "celltype",
    values_to = "af"
  ) |>
  tidyr::nest(.by = c(variant, variant_type, celltype)) |>
  dplyr::mutate(t = purrr::map(data, fn_ttest_af)) -> admeta_af_ttest

admeta_af_ttest |>
  dplyr::select(-data) |>
  tidyr::unnest(t) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    rank = plog10p * abs(estimate)
  ) |>
  dplyr::arrange(dplyr::desc(rank)) |>
  dplyr::rename(ad = "Alzheimer's Disease") |>
  dplyr::filter(ad >= 5, Healthy >= 5) -> admeta_af_ttest_sig

# Save  -------------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
