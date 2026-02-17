#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-17 15:00:00
# @DESCRIPTION: Summary plots for somatic variants prediction analysis
# @VERSION: v0.1.0

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library -----------------------------------------------------------------
suppressMessages({
  load_pkg(jutils)
  load_pkg(dplyr, data.table, ggplot2, patchwork)
})

# Args --------------------------------------------------------------------
VERSION = "v0.1.0"

GetoptLong.options(help_style = "two-column")

verbose = FALSE
GetoptLong(
  "verbose!",
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
# source("00-colors.R")

# Load data ---------------------------------------------------------------

# Set paths
basedir <- path("/liulab/chunjie/data/scMOCHA")
outdir <- path(basedir, "analysis/zzz/MANUSCRIPTFIGURES")
outdirnotuse <- path(basedir, "analysis/zzz/MANUSCRIPTFIGURES-notuse")

# Check if running from high-res directory
if (basename(getwd()) == "high-res") {
  outdir <- path("../high-res-MANUSCRIPTFIGURES")
  outdirnotuse <- path("../high-res-MANUSCRIPTFIGURES-notuse")
}

# Load predictions
somatic_with_pred <- import(outdir / "SOMATIC-VARIANTS-WITH-PREDICTIONS.xlsx")

# Clean column names for easier use
names(somatic_with_pred) <- gsub(" ", "_", tolower(names(somatic_with_pred)))

log_info("Loaded {nrow(somatic_with_pred)} variant records")

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

# 1. Prediction Coverage Plot
coverage_data <- data.table(
  tool = c("APOGEE2", "nAPOGEE", "AlphaMissense", "CADD", "PolyPhen2", "SIFT"),
  count = c(
    sum(!is.na(somatic_with_pred$apogee2_score)),
    sum(!is.na(somatic_with_pred$napogee_score)),
    sum(!is.na(somatic_with_pred$alphamissense_score)),
    sum(!is.na(somatic_with_pred$cadd_phred_score)),
    sum(!is.na(somatic_with_pred$polyphen2_score)),
    sum(!is.na(somatic_with_pred$sift_score))
  )
) |>
  mutate(
    percentage = count / nrow(distinct(somatic_with_pred, variant)) * 100,
    tool = factor(tool, levels = tool[order(count, decreasing = TRUE)])
  )

p1 <- ggplot(coverage_data, aes(x = tool, y = percentage, fill = tool)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%d\n(%.1f%%)", count, percentage)),
    vjust = -0.3,
    size = 3.5
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Prediction Tool Coverage",
    subtitle = sprintf(
      "Total: %d somatic variants",
      nrow(distinct(somatic_with_pred, variant))
    ),
    x = NULL,
    y = "Coverage (%)"
  ) +
  ylim(0, 60) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 2. APOGEE2 Classification Distribution
apogee2_class <- somatic_with_pred |>
  filter(!is.na(apogee2)) |>
  distinct(variant, .keep_all = TRUE) |>
  count(apogee2) |>
  mutate(
    apogee2 = factor(
      apogee2,
      levels = c(
        "Pathogenic",
        "Likely-pathogenic",
        "VUS+",
        "VUS",
        "VUS-",
        "Likely-benign",
        "Benign"
      )
    ),
    category = case_when(
      apogee2 %in% c("Pathogenic", "Likely-pathogenic") ~ "Pathogenic",
      apogee2 %in% c("VUS+", "VUS", "VUS-") ~ "VUS",
      TRUE ~ "Benign"
    )
  ) |>
  arrange(apogee2)

# Define colors
colors_class <- c(
  "Pathogenic" = "#d62728",
  "Likely-pathogenic" = "#ff7f0e",
  "VUS+" = "#bcbd22",
  "VUS" = "#17becf",
  "VUS-" = "#9467bd",
  "Likely-benign" = "#8c564b",
  "Benign" = "#2ca02c"
)

p2 <- ggplot(apogee2_class, aes(x = apogee2, y = n, fill = apogee2)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = colors_class) +
  labs(
    title = "APOGEE2 Classification",
    subtitle = sprintf(
      "n = %d variants with predictions",
      sum(apogee2_class$n)
    ),
    x = NULL,
    y = "Count"
  ) +
  ylim(0, max(apogee2_class$n) * 1.15) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 3. Top Pathogenic Variants by Score
top_pathogenic <- somatic_with_pred |>
  filter(apogee2 %in% c("Pathogenic", "Likely-pathogenic")) |>
  distinct(variant, .keep_all = TRUE) |>
  arrange(desc(apogee2_score)) |>
  head(12) |>
  mutate(
    variant_label = sprintf("%s (%s)", variant, gene_symbol),
    variant_label = factor(variant_label, levels = rev(variant_label))
  )

p3 <- ggplot(
  top_pathogenic,
  aes(x = variant_label, y = apogee2_score, fill = apogee2)
) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = colors_class) +
  labs(
    title = "Top Pathogenic Variants",
    subtitle = "APOGEE2 scores",
    x = NULL,
    y = "APOGEE2 Score",
    fill = "Classification"
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# 4. Gene Distribution of Pathogenic Variants
gene_dist <- somatic_with_pred |>
  filter(apogee2 %in% c("Pathogenic", "Likely-pathogenic")) |>
  distinct(variant, .keep_all = TRUE) |>
  count(gene_symbol, apogee2) |>
  arrange(desc(n)) |>
  mutate(
    gene_symbol = factor(gene_symbol, levels = unique(gene_symbol))
  )

p4 <- ggplot(gene_dist, aes(x = gene_symbol, y = n, fill = apogee2)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = colors_class) +
  labs(
    title = "Genes with Pathogenic Variants",
    subtitle = "Distribution by gene",
    x = NULL,
    y = "Number of Variants",
    fill = "Classification"
  ) +
  ylim(0, max(gene_dist$n) * 1.2) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 5. Prediction Score Comparison for Pathogenic Variants
score_comparison <- top_pathogenic |>
  select(variant, gene_symbol, apogee2_score, cadd_phred_score) |>
  mutate(
    # Normalize scores to 0-1 range for comparison
    apogee2_norm = apogee2_score,
    cadd_norm = pmin(cadd_phred_score / 30, 1), # CADD max ~30
    variant_label = sprintf("%s (%s)", variant, gene_symbol)
  ) |>
  tidyr::pivot_longer(
    cols = c(apogee2_norm, cadd_norm),
    names_to = "score_type",
    values_to = "normalized_score"
  ) |>
  mutate(
    score_type = case_when(
      score_type == "apogee2_norm" ~ "APOGEE2",
      score_type == "cadd_norm" ~ "CADD (normalized)"
    ),
    variant_label = factor(
      variant_label,
      levels = sprintf(
        "%s (%s)",
        top_pathogenic$variant,
        top_pathogenic$gene_symbol
      )
    )
  )

p5 <- ggplot(
  score_comparison,
  aes(x = variant_label, y = normalized_score, fill = score_type)
) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("APOGEE2" = "#1f77b4", "CADD (normalized)" = "#ff7f0e")
  ) +
  labs(
    title = "Multi-Tool Score Comparison",
    subtitle = "Top pathogenic variants",
    x = NULL,
    y = "Normalized Score",
    fill = "Prediction Tool"
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# 6. Score Distribution Histogram
score_dist <- somatic_with_pred |>
  distinct(variant, .keep_all = TRUE) |>
  filter(!is.na(apogee2_score))

p6 <- ggplot(score_dist, aes(x = apogee2_score, fill = apogee2)) +
  geom_histogram(bins = 30, alpha = 0.8, color = "white") +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "gray30",
    linewidth = 1
  ) +
  scale_fill_manual(values = colors_class, na.value = "gray70") +
  labs(
    title = "APOGEE2 Score Distribution",
    subtitle = "All somatic variants with predictions",
    x = "APOGEE2 Score",
    y = "Count",
    fill = "Classification"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# Combine plots -----------------------------------------------------------

# Layout: 2x3 grid
layout <- "
AABBCC
DDEEFF
"

p_combined <- p1 +
  p2 +
  p4 +
  p3 +
  p5 +
  p6 +
  plot_layout(design = layout) +
  plot_annotation(
    title = "Somatic Variants Prediction Analysis Summary",
    subtitle = sprintf(
      "Analysis of %d somatic variants with multiple prediction tools",
      nrow(distinct(somatic_with_pred, variant))
    ),
    caption = "Data sources: MitImpact DB 3.1.3 (APOGEE2), nAPOGEE v1.0.0",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save --------------------------------------------------------------------

ggsave(
  filename = outdir / "11.03-SOMATIC-VARIANTS-PREDICTION-SUMMARY.pdf",
  plot = p_combined,
  width = 16,
  height = 10,
  dpi = 300
)

log_success("Saved summary plot: 11.03-SOMATIC-VARIANTS-PREDICTION-SUMMARY.pdf")

ggsave(
  filename = outdir / "11.03-SOMATIC-VARIANTS-PREDICTION-SUMMARY.png",
  plot = p_combined,
  width = 16,
  height = 10,
  dpi = 300
)

log_success("Saved summary plot: 11.03-SOMATIC-VARIANTS-PREDICTION-SUMMARY.png")

# Individual plots for manuscript
ggsave(
  filename = outdir / "11.03a-prediction-coverage.pdf",
  plot = p1,
  width = 6,
  height = 5
)

ggsave(
  filename = outdir / "11.03b-apogee2-classification.pdf",
  plot = p2,
  width = 6,
  height = 5
)

ggsave(
  filename = outdir / "11.03c-top-pathogenic-variants.pdf",
  plot = p3,
  width = 8,
  height = 6
)

ggsave(
  filename = outdir / "11.03d-gene-distribution.pdf",
  plot = p4,
  width = 7,
  height = 5
)

log_success("Saved individual plots (a-d)")

# Session info ------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
