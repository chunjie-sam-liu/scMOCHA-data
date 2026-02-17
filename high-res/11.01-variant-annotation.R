#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-13 03:24:17
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

# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv(".env")
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})
outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
apogee2 <- import(
  outdirnotuse / "MitImpact_db_3.1.3.txt",
  format = "tsv",
  lazy = FALSE
) |>
  mutate(variant = glue("{Start}{Ref}>{Alt}"))

apogee2_n <- import(
  outdirnotuse / "nAPOGEE_v1.0.0.txt",
  format = "tsv",
  lazy = FALSE
) |>
  mutate(variant = glue("{start}{ref}>{alt}"))

ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::mutate(
    coord = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        # .v <- gse_data_variant_classification_clusteraf_bulkaf$variant[[1]]
        pos <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\1",
          x = .v
        ) |>
          as.integer()
        ref <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\2",
          x = .v
        )
        alt <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\3",
          x = .v
        )
        data.table(
          seqnames = "MT",
          start = pos,
          end = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  )

variant_annotation <- import(
  outdir / "VARIANT-ANNOTATION-TABLE.xlsx"
)

SOMATIC_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "somatic")
HOMO_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "homo")
HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "hete")
HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))
HAPLO <- ALLVARIANTS |>
  dplyr::filter(variant_type == "haplo")

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

ggvenn::ggvenn(
  list(
    "APOGEE2" = apogee2$variant,
    "nAPOGEE" = apogee2_n$variant,
    "ALLVARIANTS" = HOMO_HETE_VARIANTS$variant
  ),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5,
  set_name_size = 4
) +
  theme(
    legend.position = "none"
  ) -> p_homo_hete

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "variant_annotation_homo_hete_venn.pdf",
  plot = p_homo_hete,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

data.table(
  variant = c(HOMO_HETE_VARIANTS$variant, apogee2$variant, apogee2_n$variant) |>
    sort() |>
    unique()
) |>
  mutate(
    `HOMO+HETE` = if_else(variant %in% HOMO_HETE_VARIANTS$variant, TRUE, FALSE),
    `HOMO` = if_else(variant %in% HOMO_VARIANTS$variant, TRUE, FALSE),
    `HETE` = if_else(variant %in% HETE_VARIANTS$variant, TRUE, FALSE),
    `SOMATIC` = if_else(variant %in% SOMATIC_VARIANTS$variant, TRUE, FALSE),
    `HAPLO` = if_else(variant %in% HAPLO$variant, TRUE, FALSE),
    APOGEE2 = if_else(variant %in% apogee2$variant, TRUE, FALSE),
    APOGEE2_N = if_else(variant %in% apogee2_n$variant, TRUE, FALSE)
  ) -> foreuler

library(eulerr)
color_variant_type <- c(
  "HOMO+HETE" = "#3AF712",
  "APOGEE2" = "#0073C2FF",
  "APOGEE2_N" = "#EFC000FF",
  "HOMO" = "#8DD3C7FF",
  "HETE" = "#FFFFB3FF",
  "SOMATIC" = "red",
  "HAPLO" = "#BEBADAFF"
)


{
  plot(
    euler(
      foreuler[, c(
        "HOMO+HETE",
        "APOGEE2",
        "APOGEE2_N",
        # "HOMO",
        # "HETE",
        # "SOMATIC",
        "HAPLO"
      )],
      # shape = "ellipse",
      control = list(extraopt = FALSE)
    ),
    # quantities = list(type = c("counts", "percent"), font = 5),
    quantities = list(type = c("counts"), font = 5),
    labels = list(fontfamily = "serif"),
    edges = list(lty = 3),
    # fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
    fills = color_variant_type
  ) -> p_euler_homo_hete
  pdf(
    path(
      outdirnotuse / "variant-annotation",
      "Variant-type-Euler-HOMOHETE.pdf"
    ),
    width = 8,
    height = 6
  )
  print(p_euler_homo_hete)
  dev.off()
}


{
  plot(
    euler(
      foreuler[, c(
        # "HOMO+HETE",
        "APOGEE2",
        "APOGEE2_N",
        # "HOMO",
        # "HETE",
        "SOMATIC"
        # "HAPLO"
      )],
      # shape = "ellipse",
      control = list(extraopt = FALSE)
    ),
    # quantities = list(type = c("counts", "percent"), font = 5),
    quantities = list(type = c("counts"), font = 5),
    labels = list(fontfamily = "serif"),
    edges = list(lty = 3),
    # fills = c("#BEBADAFF", "#8DD3C7FF", "#FFFFB3FF", "red")
    fills = color_variant_type
  ) -> p_euler_somatic
  pdf(
    path(
      outdirnotuse / "variant-annotation",
      "Variant-type-Euler-SOMATIC.pdf"
    ),
    width = 8,
    height = 6
  )
  print(p_euler_somatic)
  dev.off()
}


#
#
# ? variant annotation --------------------------------------------------------------------
#
#

variant_annotation |> count(Locus) |> arrange(-n) |> print(n = Inf)
variant_annotation |> count(aachange) |> arrange(-n)
variant_annotation |>
  mutate(
    aachange_new = if_else(
      is.na(aachange),
      "Unknown",
      aachange
    )
  ) |>
  mutate(
    aachange_new = if_else(
      grepl("MitoTIP", aachange_new),
      "MitoTIP",
      aachange_new
    )
  ) |>
  mutate(
    Locus_new = if_else(
      is.na(Locus),
      "Unknown",
      Locus
    )
  ) |>
  mutate(
    # Check if variant is in control region/D-loop
    is_dloop = grepl(
      "ATT,CR:|CR:Control Region|CR:HVS|CR:OH|CR:7S|CR:CSB|CR:TF|CR:mt|CR:HPR|CR:PL|CR:PH",
      Locus_new
    )
  ) |>
  mutate(
    aachange_group = case_when(
      grepl(":", aachange_new) ~ "Protein-coding",
      grepl("rRNA", aachange_new) ~ "rRNA",
      # MitoTIP scores in D-loop should be classified as non-coding, not tRNA
      aachange_new == "MitoTIP" & is_dloop ~ "Non-coding",
      aachange_new == "MitoTIP" ~ "tRNA",
      aachange_new == "Unknown" ~ "Unknown",
      aachange_new == "non-coding" ~ "Non-coding"
    )
  ) |>
  mutate(
    Locus_group = case_when(
      # Protein-coding genes
      grepl("^ND[1-6]$|^ND4L$|^ND3$", Locus_new) ~ "ND genes",
      grepl(
        "^Cytb$|^COIII$|^COI$|^COII$",
        Locus_new
      ) ~ "CO/Cytb genes",
      grepl("^ATPase[68]$", Locus_new) ~ "ATPase genes",
      grepl("^ND\\d,ND\\d", Locus_new) ~ "ND genes",
      grepl("^ATPase6,ATPase8$", Locus_new) ~ "ATPase genes",
      # rRNA genes
      grepl("^12S$|^16S$|^12S,|^16S,|-,12S|-,16S", Locus_new) ~ "rRNA",
      # tRNA genes (single letters or L()/S() notation)
      grepl(
        "^[A-Z]$|^[A-Z],|,[A-Z]$|^L\\(|^S\\(|^-,[A-Z]|^OL,[A-Z]",
        Locus_new
      ) ~ "tRNA",
      # Control region (all ATT,CR: and CR: entries)
      grepl(
        "ATT,CR:|CR:Control Region|CR:HVS|CR:OH|CR:7S|CR:CSB|CR:TF|CR:mt|CR:HPR|CR:PL|CR:PH",
        Locus_new
      ) ~ "D-loop",
      # Non-coding regions
      grepl("^NC[0-9]+$|^ATT,NC", Locus_new) ~ "Non-coding (intergenic)",
      # Origin of replication
      grepl("^OL$", Locus_new) ~ "Origin of Replication",
      # Unknown
      Locus_new == "Unknown" ~ "Unknown",
      # Everything else
      TRUE ~ "Other"
    )
  ) |>
  select(-is_dloop) |>
  # Add APOGEE2 pathogenicity predictions for protein-coding variants
  left_join(
    apogee2 |>
      select(
        variant,
        APOGEE2_score,
        APOGEE2_probability,
        APOGEE2
      ),
    by = "variant"
  ) |>
  # Add nAPOGEE pathogenicity predictions for non-coding variants
  left_join(
    apogee2_n |>
      select(
        variant,
        nAPOGEE_score,
        nAPOGEE_posterior_probability,
        pathogenicity_assessment
      ),
    by = "variant",
    suffix = c("", "_napogee")
  ) -> variant_annotation_grouped

export(
  variant_annotation_grouped |>
    mutate(
      prediction_class = case_when(
        aachange_group == "Protein-coding" & !is.na(APOGEE2) ~ APOGEE2,
        aachange_group == "Protein-coding" & is.na(APOGEE2) ~ "No prediction",
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          !is.na(pathogenicity_assessment) ~ pathogenicity_assessment,
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          is.na(pathogenicity_assessment) ~ "No prediction",
        TRUE ~ "Not applicable"
      ),
      # Standardize classification names
      prediction_class = case_when(
        prediction_class %in% c("pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c("likely pathogenic", "Likely-pathogenic") ~ "Likely pathogenic",
        prediction_class %in% c("benign", "Benign") ~ "Benign",
        prediction_class %in%
          c("likely benign", "Likely-benign") ~ "Likely benign",
        TRUE ~ prediction_class
      ),
      prediction_category = case_when(
        aachange_group == "Protein-coding" ~ "APOGEE2",
        aachange_group %in% c("Non-coding", "tRNA", "rRNA") ~ "nAPOGEE",
        TRUE ~ "Not applicable"
      )
    ),
  outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

export(
  variant_annotation_grouped,
  outdirnotuse / "variant-annotation" / "VARIANT-ANNOTATION-TABLE-APOGEE2.qs",
)

# Check the grouping
variant_annotation_grouped |>
  count(Locus_group) |>
  arrange(-n) |>
  print()

variant_annotation_grouped |>
  count(aachange_group) |>
  arrange(-n) |>
  print()

# Check APOGEE2 and nAPOGEE coverage
log_info("APOGEE2 coverage for protein-coding variants...")
variant_annotation_grouped |>
  filter(aachange_group == "Protein-coding") |>
  summarise(
    total = n(),
    with_APOGEE2 = sum(!is.na(APOGEE2)),
    pct = with_APOGEE2 / total * 100
  ) |>
  print()

log_info("nAPOGEE coverage for non-coding variants...")
variant_annotation_grouped |>
  filter(aachange_group == "Non-coding") |>
  summarise(
    total = n(),
    with_nAPOGEE = sum(!is.na(pathogenicity_assessment)),
    pct = with_nAPOGEE / total * 100
  ) |>
  print()

log_info("APOGEE2 classifications...")
variant_annotation_grouped |>
  count(APOGEE2) |>
  arrange(-n) |>
  print()

log_info("APOGEE2 classifications for protein-coding variants specifically...")
variant_annotation_grouped |>
  filter(aachange_group == "Protein-coding") |>
  count(APOGEE2) |>
  arrange(-n) |>
  print()

log_info("nAPOGEE classifications...")
variant_annotation_grouped |>
  count(pathogenicity_assessment) |>
  arrange(-n) |>
  print()

# Plot Locus_group vs aachange_group with counts
load_pkg(ggplot2)

variant_annotation_grouped |>
  count(Locus_group, aachange_group) |>
  mutate(
    total = sum(n),
    pct = n / total * 100,
    label = sprintf("%d\n(%.1f%%)", n, pct),
    Locus_group = factor(
      Locus_group,
      levels = c(
        "ND genes",
        "CO/Cytb genes",
        "ATPase genes",
        "rRNA",
        "tRNA",
        "D-loop",
        "Non-coding (intergenic)",
        "Origin of Replication",
        "Unknown",
        "Other"
      )
    ),
    aachange_group = factor(
      aachange_group,
      levels = c(
        "Protein-coding",
        "rRNA",
        "tRNA",
        "Non-coding",
        "Unknown"
      )
    )
  ) |>
  ggplot(aes(x = aachange_group, y = Locus_group, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.5, color = "black") +
  scale_fill_gradient(
    low = "#f7f7f7",
    high = "#d62728",
    name = "Count"
  ) +
  labs(
    x = "Amino Acid Change Group",
    y = "Locus Group",
    title = "Variant Distribution by Locus and AA Change Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) -> p_locus_aachange

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "locus_group_vs_aachange_group.pdf",
  plot = p_locus_aachange,
  width = 10,
  height = 8,
  dpi = 300
)

# Plot APOGEE prediction coverage ----------------------------------------

#' Create stacked bar chart for APOGEE pathogenicity predictions
#'
#' @param data Variant annotation data with aachange_group and APOGEE predictions
#' @param title Plot title
#' @param verbose Print plot data summary
#'
#' @return ggplot object
create_apogee_coverage_plot <- function(
  data,
  title = "APOGEE Pathogenicity Predictions",
  verbose = FALSE
) {
  # Define colors
  color_prediction <- c(
    "Pathogenic" = "#d62728",
    "Likely pathogenic" = "#fc7373",
    "VUS+" = "#ffbb78",
    "VUS" = "#ffed6f",
    "VUS-" = "#fff7bc",
    "Likely benign" = "#98df8a",
    "Benign" = "#2ca02c",
    "No prediction" = "#d3d3d3"
  )

  # Prepare plot data
  plot_data <- data |>
    mutate(
      prediction_class = case_when(
        aachange_group == "Protein-coding" & !is.na(APOGEE2) ~ APOGEE2,
        aachange_group == "Protein-coding" & is.na(APOGEE2) ~ "No prediction",
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          !is.na(pathogenicity_assessment) ~ pathogenicity_assessment,
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          is.na(pathogenicity_assessment) ~ "No prediction",
        TRUE ~ "Not applicable"
      ),
      # Standardize classification names
      prediction_class = case_when(
        prediction_class %in% c("pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c("likely pathogenic", "Likely-pathogenic") ~ "Likely pathogenic",
        prediction_class %in% c("benign", "Benign") ~ "Benign",
        prediction_class %in%
          c("likely benign", "Likely-benign") ~ "Likely benign",
        TRUE ~ prediction_class
      ),
      prediction_category = case_when(
        aachange_group == "Protein-coding" ~ "APOGEE2",
        aachange_group %in% c("Non-coding", "tRNA", "rRNA") ~ "nAPOGEE",
        TRUE ~ "Not applicable"
      )
    ) |>
    filter(prediction_category != "Not applicable") |>
    count(aachange_group, prediction_class) |>
    group_by(aachange_group) |>
    mutate(
      total = sum(n),
      pct = n / total * 100,
      label = as.character(n)
    ) |>
    ungroup() |>
    mutate(
      aachange_group = factor(
        aachange_group,
        levels = c("Protein-coding", "rRNA", "tRNA", "Non-coding")
      ),
      prediction_class = factor(
        prediction_class,
        levels = c(
          "Pathogenic",
          "Likely pathogenic",
          "VUS+",
          "VUS",
          "VUS-",
          "Likely benign",
          "Benign",
          "No prediction"
        )
      )
    )

  # Print summary if verbose
  if (verbose) {
    log_info("Plot data summary:")
    print(plot_data, n = Inf)
  }

  # Create plot
  p <- plot_data |>
    ggplot(aes(x = aachange_group, y = n, fill = prediction_class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(
      values = color_prediction,
      name = "Pathogenicity"
    ) +
    labs(
      x = "Variant Type",
      y = "Number of Variants",
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  return(p)
}

# Generate APOGEE coverage plot
p_apogee_coverage <- create_apogee_coverage_plot(
  data = variant_annotation_grouped,
  title = "APOGEE Pathogenicity Predictions",
  verbose = TRUE
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "apogee_pathogenicity_predictions.pdf",
  plot = p_apogee_coverage,
  width = 10,
  height = 6,
  dpi = 300
)

# Two-layer pie chart: aachange_group (inner) + APOGEE prediction (outer) ----

#' Create two-layer pie chart for variant type and pathogenicity
#'
#' @param data Variant annotation data with aachange_group and prediction columns
#' @param title Plot title
#'
#' @return ggplot object
create_twolayer_pie <- function(
  data,
  title = "Two-Layer Distribution: Variant Type and Pathogenicity"
) {
  # Load required libraries
  library(ggnewscale)
  library(ggrepel)

  # Define colors - distinct from APOGEE colors
  color_aachange = c(
    "Protein-coding" = "#7570b3",
    "rRNA" = "#4292c6",
    "tRNA" = "#41b6c4",
    "Non-coding" = "#e7a77c",
    "Unknown" = "#bdbdbd"
  )

  color_prediction = c(
    "Pathogenic" = "#d62728",
    "Likely pathogenic" = "#fc7373",
    "VUS+" = "#ffbb78",
    "VUS" = "#ffed6f",
    "VUS-" = "#fff7bc",
    "Likely benign" = "#98df8a",
    "Benign" = "#2ca02c",
    "No prediction" = "#d3d3d3"
  )

  # Prepare pie data
  pie_data <- data |>
    mutate(
      prediction_class = case_when(
        aachange_group == "Protein-coding" & !is.na(APOGEE2) ~ APOGEE2,
        aachange_group == "Protein-coding" & is.na(APOGEE2) ~ "No prediction",
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          !is.na(pathogenicity_assessment) ~ pathogenicity_assessment,
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          is.na(pathogenicity_assessment) ~ "No prediction",
        TRUE ~ "Not applicable"
      ),
      # Standardize classification names
      prediction_class = case_when(
        prediction_class %in% c("pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c("likely pathogenic", "Likely-pathogenic") ~ "Likely pathogenic",
        prediction_class %in% c("benign", "Benign") ~ "Benign",
        prediction_class %in%
          c("likely benign", "Likely-benign") ~ "Likely benign",
        TRUE ~ prediction_class
      )
    ) |>
    filter(prediction_class != "Not applicable") |>
    count(aachange_group, prediction_class) |>
    mutate(
      aachange_group = factor(
        aachange_group,
        levels = c("Protein-coding", "rRNA", "tRNA", "Non-coding", "Unknown")
      ),
      prediction_class = factor(
        prediction_class,
        levels = c(
          "Pathogenic",
          "Likely pathogenic",
          "VUS+",
          "VUS",
          "VUS-",
          "Likely benign",
          "Benign",
          "No prediction"
        )
      )
    )

  # Calculate positions for outer ring
  pie_data_positioned <- pie_data |>
    group_by(aachange_group) |>
    mutate(
      group_total = sum(n),
      group_pct = n / group_total * 100
    ) |>
    ungroup() |>
    mutate(
      total = sum(n),
      overall_pct = n / total * 100
    )

  # Create inner layer data (aachange_group totals)
  inner_data <- pie_data_positioned |>
    group_by(aachange_group) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(
      fraction = n / sum(n),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      percentage = sprintf("%.1f%%", fraction * 100),
      label = paste0(as.character(aachange_group), "\n", percentage)
    )

  # Create outer layer data (prediction_class by aachange_group)
  outer_data <- pie_data_positioned |>
    arrange(aachange_group, prediction_class) |>
    group_by(aachange_group) |>
    mutate(
      group_fraction = n / sum(n),
      group_ymax = cumsum(group_fraction),
      group_ymin = c(0, head(group_ymax, n = -1))
    ) |>
    ungroup() |>
    left_join(
      inner_data |>
        select(aachange_group, ymin_inner = ymin, ymax_inner = ymax),
      by = "aachange_group"
    ) |>
    mutate(
      ymin = ymin_inner + (ymax_inner - ymin_inner) * group_ymin,
      ymax = ymin_inner + (ymax_inner - ymin_inner) * group_ymax,
      labelPosition = (ymax + ymin) / 2,
      percentage = sprintf("%.1f%%", group_fraction * 100),
      label = paste0(as.character(prediction_class), "\n", percentage)
    )

  # Create the plot
  p <- ggplot() +
    # Inner layer (aachange_group)
    geom_rect(
      data = inner_data,
      aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = aachange_group),
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      data = inner_data,
      aes(x = 3.5, y = labelPosition, label = label),
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = color_aachange,
      name = "Variant Type",
      guide = guide_legend(order = 1)
    ) +
    new_scale_fill() +
    # Outer layer (prediction_class)
    geom_rect(
      data = outer_data,
      aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 5,
        xmin = 4,
        fill = prediction_class
      ),
      color = "white",
      linewidth = 0.3
    ) +
    geom_text_repel(
      data = outer_data,
      aes(x = 4.5, y = labelPosition, label = label),
      size = 3,
      lineheight = 0.85,
      segment.size = 0.3,
      segment.color = "grey40",
      min.segment.length = 0,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 3,
      max.overlaps = Inf,
      direction = "both",
      seed = 42
    ) +
    scale_fill_manual(
      values = color_prediction,
      name = "Pathogenicity",
      guide = guide_legend(order = 2)
    ) +
    coord_polar(theta = "y") +
    xlim(c(2, 5)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    ) +
    labs(title = title)

  return(p)
}

#' Create one-layer pie chart for pathogenicity predictions
#'
#' @param data Variant annotation data with prediction columns
#' @param title Plot title
#'
#' @return ggplot object
create_onelayer_pie <- function(
  data,
  title = "Pathogenicity Distribution"
) {
  # Simplified color palette
  color_prediction_simple <- c(
    "Pathogenic" = "#ff0000",
    "Likely pathogenic" = "#fc7373",
    "VUS" = "#ffed6f",
    "Likely benign" = "#98df8a",
    "Benign" = "#2ca02c",
    "No prediction" = "#d3d3d3"
  )

  # Prepare pie data with grouped categories
  pie_data <- data |>
    mutate(
      prediction_class = case_when(
        aachange_group == "Protein-coding" & !is.na(APOGEE2) ~ APOGEE2,
        aachange_group == "Protein-coding" & is.na(APOGEE2) ~ "No prediction",
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          !is.na(pathogenicity_assessment) ~ pathogenicity_assessment,
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          is.na(pathogenicity_assessment) ~ "No prediction",
        TRUE ~ "Not applicable"
      ),
      # Standardize and group classifications
      prediction_class_grouped = case_when(
        prediction_class %in% c("pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c("likely pathogenic", "Likely-pathogenic") ~ "Likely pathogenic",
        prediction_class %in% c("VUS-", "VUS", "VUS+") ~ "VUS",
        prediction_class %in% c("benign", "Benign") ~ "Benign",
        prediction_class %in%
          c("likely benign", "Likely-benign") ~ "Likely benign",
        prediction_class == "No prediction" ~ "No prediction",
        TRUE ~ "Not applicable"
      )
    ) |>
    filter(prediction_class_grouped != "Not applicable") |>
    count(prediction_class_grouped) |>
    mutate(
      prediction_class_grouped = factor(
        prediction_class_grouped,
        levels = c(
          "Pathogenic",
          "Likely pathogenic",
          "VUS",
          "Likely benign",
          "Benign",
          "No prediction"
        )
      )
    ) |>
    arrange(prediction_class_grouped) |>
    mutate(
      fraction = n / sum(n),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      percentage = sprintf("%.1f%%", fraction * 100),
      label = paste0(prediction_class_grouped, "\n", percentage)
    )

  # Create the plot
  p <- ggplot(pie_data) +
    geom_rect(
      aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = prediction_class_grouped
      ),
      color = "white",
      linewidth = 1
    ) +
    geom_text(
      aes(x = 3.5, y = labelPosition, label = label),
      size = 4,
      fontface = "bold",
      lineheight = 0.9
    ) +
    scale_fill_manual(
      values = color_prediction_simple,
      name = "Pathogenicity",
      breaks = c(
        "Pathogenic",
        "Likely pathogenic",
        "VUS",
        "Likely benign",
        "Benign",
        "No prediction"
      )
    ) +
    coord_polar(theta = "y") +
    # xlim(c(0, 4)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(title = title)

  return(p)
}

# Generate the two-layer pie chart
p_twolayer_pie <- create_twolayer_pie(
  data = variant_annotation_grouped,
  title = "Two-Layer Distribution: Variant Type and Pathogenicity"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "two_layer_pie_aachange_apogee.pdf",
  plot = p_twolayer_pie,
  width = 12,
  height = 8,
  dpi = 300
)

# Generate the one-layer pie chart
p_onelayer_pie <- create_onelayer_pie(
  data = variant_annotation_grouped,
  title = "Pathogenicity Distribution"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "one_layer_pie_pathogenicity.pdf",
  plot = p_onelayer_pie,
  width = 10,
  height = 8,
  dpi = 300
)

# Generate plots for somatic variants ----------------------------------------
log_info("Creating plots for somatic variants...")

variant_annotation_somatic <- variant_annotation_grouped |>
  filter(variant %in% SOMATIC_VARIANTS$variant)

log_info("Found {nrow(variant_annotation_somatic)} somatic variants")

# Two-layer pie chart for somatic variants
p_twolayer_pie_somatic <- create_twolayer_pie(
  data = variant_annotation_somatic,
  title = "Somatic Variants: Type and Pathogenicity"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "two_layer_pie_aachange_apogee_somatic.pdf",
  plot = p_twolayer_pie_somatic,
  width = 12,
  height = 8,
  dpi = 300
)

# One-layer pie chart for somatic variants
p_onelayer_pie_somatic <- create_onelayer_pie(
  data = variant_annotation_somatic,
  title = "Somatic Variants: Pathogenicity Distribution"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "one_layer_pie_pathogenicity_somatic.pdf",
  plot = p_onelayer_pie_somatic,
  width = 10,
  height = 8,
  dpi = 300
)

# APOGEE coverage plot for somatic variants
p_apogee_coverage_somatic <- create_apogee_coverage_plot(
  data = variant_annotation_somatic,
  title = "APOGEE Pathogenicity Predictions Coverage for Somatic Variants"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "apogee_pathogenicity_predictions_somatic.pdf",
  plot = p_apogee_coverage_somatic,
  width = 10,
  height = 6,
  dpi = 300
)

# Generate plots for heteroplasmic variants ----------------------------------
log_info("Creating plots for heteroplasmic variants...")

variant_annotation_hete <- variant_annotation_grouped |>
  filter(variant %in% HETE_VARIANTS$variant)

log_info("Found {nrow(variant_annotation_hete)} heteroplasmic variants")

# Two-layer pie chart for heteroplasmic variants
p_twolayer_pie_hete <- create_twolayer_pie(
  data = variant_annotation_hete,
  title = "Heteroplasmic Variants: Type and Pathogenicity"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "two_layer_pie_aachange_apogee_hete.pdf",
  plot = p_twolayer_pie_hete,
  width = 12,
  height = 8,
  dpi = 300
)

# One-layer pie chart for heteroplasmic variants
p_onelayer_pie_hete <- create_onelayer_pie(
  data = variant_annotation_hete,
  title = "Heteroplasmic Variants: Pathogenicity Distribution"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "one_layer_pie_pathogenicity_hete.pdf",
  plot = p_onelayer_pie_hete,
  width = 10,
  height = 8,
  dpi = 300
)

# APOGEE coverage plot for heteroplasmic variants
p_apogee_coverage_hete <- create_apogee_coverage_plot(
  data = variant_annotation_hete,
  title = "APOGEE Pathogenicity Predictions Coverage for Heteroplasmic Variants"
)

ggsave(
  filename = outdirnotuse /
    "variant-annotation" /
    "apogee_pathogenicity_predictions_hete.pdf",
  plot = p_apogee_coverage_hete,
  width = 10,
  height = 6,
  dpi = 300
)


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
