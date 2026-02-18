#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-18 15:47:28
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

variant_annotation_grouped <- import(
  outdirnotuse / "variant-annotation" / "VARIANT-ANNOTATION-TABLE-APOGEE2.qs",
)
# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------
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
        # Protein-coding: APOGEE2 prediction available
        aachange_group == "Protein-coding" & !is.na(APOGEE2) ~ APOGEE2,
        # Protein-coding: stop-gain (TAA/TGA stop codon) without APOGEE2 → Likely pathogenic
        aachange_group == "Protein-coding" &
          is.na(APOGEE2) &
          grepl("\\*", aachange) ~ "Likely pathogenic",
        # Protein-coding: synonymous without APOGEE2 → Likely benign
        aachange_group == "Protein-coding" &
          is.na(APOGEE2) &
          grepl("^\\w+:[A-Z]\\d+[A-Z]$", aachange) &
          substr(aachange, nchar(aachange), nchar(aachange)) ==
            substr(
              aachange,
              regexpr(":", aachange) + 1,
              regexpr(":", aachange) + 1
            ) ~
          "Likely benign",
        # Protein-coding: other without APOGEE2 → No prediction
        aachange_group == "Protein-coding" & is.na(APOGEE2) ~ "No prediction",
        # Non-coding/tRNA/rRNA: pathogenicity assessment available
        aachange_group %in%
          c("Non-coding", "tRNA", "rRNA") &
          !is.na(pathogenicity_assessment) ~ pathogenicity_assessment,
        # Non-coding without assessment → VUS
        aachange_group == "Non-coding" &
          is.na(pathogenicity_assessment) ~ "VUS",
        # tRNA/rRNA without assessment → No prediction
        aachange_group %in%
          c("tRNA", "rRNA") &
          is.na(pathogenicity_assessment) ~ "No prediction",
        TRUE ~ "Not applicable"
      ),
      # Standardize and group classifications
      prediction_class_grouped = case_when(
        prediction_class %in% c("pathogenic", "Pathogenic") ~ "Pathogenic",
        prediction_class %in%
          c(
            "likely pathogenic",
            "Likely-pathogenic",
            "Likely pathogenic"
          ) ~ "Likely pathogenic",
        prediction_class %in% c("VUS-", "VUS", "VUS+") ~ "VUS",
        prediction_class %in% c("benign", "Benign") ~ "Benign",
        prediction_class %in%
          c(
            "likely benign",
            "Likely-benign",
            "Likely benign"
          ) ~ "Likely benign",
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

# Main --------------------------------------------------------------------

# Generate plots for somatic variants ----------------------------------------
log_info("Creating plots for somatic variants...")

variant_annotation_somatic <- variant_annotation_grouped |>
  filter(variant %in% SOMATIC_VARIANTS$variant)


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


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
