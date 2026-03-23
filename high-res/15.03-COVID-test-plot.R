#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: Cluster-level COVID-19 variant t-test scatter and violin plots

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args ---------------------------------------------------------------------

VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)

# Logger -------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}

# Load data ----------------------------------------------------------------

load_pkg(jutils)
dotenv()
source(path(Sys.getenv("HIGHRESDIR"), "00-colors.R"))

variant_annotation <- import(
  path(Sys.getenv("OUTDIR")) / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"

covid_variant_ttest_cluster <- import(
  covid_outdir / "COVID19-variant-af-ttest-cluster.qs2"
)

# Packages -----------------------------------------------------------------

suppressMessages({
  load_pkg(ggh4x, ggbeeswarm, ggnewscale, ggrepel, ggsignif)
})

# Constants ----------------------------------------------------------------

color_celltype_bulk <- c("Pseudo-bulk" = "red", color_celltype)

# Functions ----------------------------------------------------------------

theme_covid_panel <- function() {
  theme(
    plot.margin = margin(t = 0.2, b = 0.1, l = 0.1, r = 0.2, unit = "cm"),
    panel.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    legend.position = "top",
    legend.key = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, color = "black", face = "bold")
  )
}

strip_celltype <- function() {
  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = color_celltype_bulk, color = NA),
    text_x = ggh4x::elem_list_text(colour = "white", face = "bold")
  )
}

normalize_celltype <- function(x) {
  x <- gsub(pattern = "_", replacement = " ", x = x)
  ifelse(tolower(x) == "bulk", "Pseudo-bulk", x)
}

rank_variants <- function(ttest_data, top_n = 8) {
  ttest_data |>
    select(-data) |>
    unnest(t) |>
    filter(!is.na(p.value)) |>
    filter(Healthy >= 5, `COVID-19` >= 5) |>
    mutate(
      rank_score = -log10(pmax(p.value, 1e-300)) * abs(estimate)
    ) -> ranked

  ranked |>
    filter(p.value < 0.05) |>
    group_by(variant) |>
    summarise(rank_score = max(rank_score, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(rank_score)) -> sig_ranked

  if (nrow(sig_ranked) == 0) {
    log_warn("No significant variants found; falling back to strongest overall hits")
    ranked |>
      group_by(variant) |>
      summarise(rank_score = max(rank_score, na.rm = TRUE), .groups = "drop") |>
      arrange(desc(rank_score)) -> sig_ranked
  }

  head(sig_ranked$variant, top_n)
}

plot_ttest_scatter <- function(ttest_data, label = "cluster") {
  ttest_data |>
    select(-data) |>
    unnest(t) |>
    filter(!is.na(p.value)) |>
    filter(Healthy >= 5, `COVID-19` >= 5) |>
    mutate(
      log10p = -log10(pmax(p.value, 1e-300)),
      celltype = normalize_celltype(celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk)),
      effect = estimate2 - estimate1,
      point_color = ifelse(p.value < 0.05 & abs(effect) > 0.03, "red", "black"),
      label = ifelse(point_color == "red", variant, NA)
    ) -> forplot

  ggplot(forplot, aes(x = effect, y = log10p)) +
    geom_point(aes(color = point_color)) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      data = forplot |> filter(!is.na(label)),
      aes(label = label),
      size = 3,
      show.legend = FALSE,
      segment.color = "black",
      max.overlaps = Inf
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = 20, color = "red") +
    geom_vline(xintercept = 0, linetype = 20, color = "red") +
    ggh4x::facet_grid2(~celltype, strip = strip_celltype()) +
    theme_covid_panel() +
    theme(axis.ticks.x = element_blank()) +
    labs(
      title = paste("COVID-19 variant t-test -", label, "level"),
      y = "-log10(p-value)",
      x = "Effect Size (COVID-19 - Healthy)"
    )
}

plot_af_violin <- function(ttest_data, variants, label = "cluster") {
  ttest_data |>
    select(variant, celltype, data) |>
    filter(variant %in% variants) |>
    unnest(data) |>
    mutate(
      celltype = normalize_celltype(celltype),
      celltype = factor(celltype, levels = names(color_celltype_bulk))
    ) |>
    ggplot(aes(x = disease)) +
    ggh4x::facet_grid2(
      variant ~ celltype,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(fill = color_celltype_bulk, color = NA),
        text_x = ggh4x::elem_list_text(colour = "white", face = "bold"),
        background_y = ggh4x::elem_list_rect(fill = "black", color = NA),
        text_y = ggh4x::elem_list_text(colour = "white", face = "bold")
      )
    ) +
    geom_violin(aes(y = af, fill = disease), alpha = 0.7, color = NA) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = color_disease, name = "Disease") +
    ggbeeswarm::geom_quasirandom(
      aes(y = af, color = disease),
      size = 1,
      dodge.width = 0.75,
      alpha = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_disease, name = "Disease") +
    ggsignif::geom_signif(
      aes(y = af),
      comparisons = list(c("COVID-19", "Healthy")),
      y_position = 0.8
    ) +
    theme_covid_panel() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = paste("COVID-19 variant AF distribution -", label, "level"),
      y = "Allele Frequency"
    )
}

# Main ---------------------------------------------------------------------

topvariants <- rank_variants(covid_variant_ttest_cluster)
log_info("Selected {length(topvariants)} top variants for COVID-19 violin plots")

topvariants_qs2 <- as.character(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2")
topvariants_qs <- as.character(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs")

if (file.exists(topvariants_qs2)) {
  invisible(file.remove(topvariants_qs2))
}

if (file.exists(topvariants_qs)) {
  invisible(file.remove(topvariants_qs))
}

export(
  topvariants,
  covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2"
)

if (file.exists(topvariants_qs)) {
  file.copy(
    topvariants_qs,
    topvariants_qs2,
    overwrite = TRUE
  )
  invisible(file.remove(topvariants_qs))
}

log_info("Plotting cluster-level COVID-19 t-test scatter")
p_scatter_cluster <- plot_ttest_scatter(
  ttest_data = covid_variant_ttest_cluster,
  label = "cluster"
)
ggsave(
  filename = covid_outdir / "COVID19-variant-af-ttest-cluster.pdf",
  plot = p_scatter_cluster,
  width = 13,
  height = 3.5
)

log_info("Plotting cluster-level COVID-19 AF violin")
p_violin_cluster <- plot_af_violin(
  ttest_data = covid_variant_ttest_cluster,
  variants = topvariants,
  label = "cluster"
)
ggsave(
  filename = covid_outdir / "COVID19-variant-af-violin-cluster.pdf",
  plot = p_violin_cluster,
  width = 16,
  height = max(8, length(topvariants) * 1.8),
  device = cairo_pdf
)

if (isTRUE(verbose)) {
  sessionInfo()
}
