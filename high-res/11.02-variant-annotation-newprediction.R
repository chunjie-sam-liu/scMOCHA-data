#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-17 14:06:03
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
source("00-colors.R")

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

variant_annotation_grouped <- import(
  outdirnotuse / "variant-annotation" / "VARIANT-ANNOTATION-TABLE-APOGEE2.qs"
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
variant_annotation_somatic <- variant_annotation_grouped |>
  filter(variant %in% SOMATIC_VARIANTS$variant) |>
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
  )

variant_annotation_somatic |>
  filter(prediction_class_grouped == "No prediction") -> no_prediction_somatic

# Explore alternative predictions for "No prediction" variants ---------
log_info(
  "Exploring {nrow(no_prediction_somatic)} variants with 'No prediction'"
)

# Merge with apogee2 to get all prediction columns
no_prediction_with_apogee2 <- no_prediction_somatic |>
  dplyr::left_join(
    apogee2 |>
      dplyr::select(
        variant,
        aachange_group = Molecule_type,
        # Prediction tools
        PolyPhen2,
        PolyPhen2_score,
        SIFT,
        SIFT_score,
        SIFT4G,
        SIFT4G_score,
        VEST,
        VEST_pvalue,
        Mitoclass1,
        SNPDryad,
        SNPDryad_score,
        MutationTaster,
        MutationTaster_score,
        FATHMM,
        FATHMM_score,
        AlphaMissense,
        AlphaMissense_score,
        CADD,
        CADD_phred_score,
        PROVEAN,
        PROVEAN_score,
        MutationAssessor,
        MutationAssessor_score,
        EFIN_SP,
        EFIN_SP_score,
        EFIN_HD,
        EFIN_HD_score,
        MLC,
        MLC_score,
        APOGEE1,
        APOGEE1_score,
        APOGEE2_alt = APOGEE2,
        APOGEE2_score,
        CAROL,
        CAROL_score,
        Condel,
        Condel_score,
        COVEC_WMV,
        COVEC_WMV_score,
        MtoolBox,
        MtoolBox_DS,
        DEOGEN2,
        DEOGEN2_score
      ),
    by = "variant",
    suffix = c("", "_apogee2")
  )

# Count available predictions per variant
prediction_tools <- c(
  "PolyPhen2",
  "SIFT",
  "SIFT4G",
  "VEST",
  "Mitoclass1",
  "SNPDryad",
  "MutationTaster",
  "FATHMM",
  "AlphaMissense",
  "CADD",
  "PROVEAN",
  "MutationAssessor",
  "EFIN_SP",
  "EFIN_HD",
  "MLC",
  "APOGEE1",
  "APOGEE2_alt",
  "CAROL",
  "Condel",
  "COVEC_WMV",
  "MtoolBox",
  "DEOGEN2"
)

no_prediction_summary <- no_prediction_with_apogee2 |>
  dplyr::mutate(
    n_predictions = rowSums(
      !is.na(across(all_of(prediction_tools))) &
        across(all_of(prediction_tools)) != "."
    )
  ) |>
  dplyr::arrange(desc(n_predictions))

log_info("Summary of alternative predictions:")
no_prediction_summary |>
  dplyr::count(aachange_group, n_predictions) |>
  dplyr::arrange(aachange_group, desc(n_predictions)) |>
  print()

# Show examples with most predictions
log_info("Examples with most alternative predictions:")
no_prediction_summary |>
  dplyr::filter(n_predictions > 0) |>
  dplyr::select(
    variant,
    aachange_group,
    n_predictions,
    PolyPhen2,
    SIFT,
    AlphaMissense,
    CADD,
    MLC,
    APOGEE2_alt
  ) |>
  head(20) |>
  print()

# Count how many prediction tools have data for each tool
log_info("Coverage of each prediction tool:")
prediction_coverage <- sapply(prediction_tools, function(tool) {
  sum(
    !is.na(no_prediction_with_apogee2[[tool]]) &
      no_prediction_with_apogee2[[tool]] != "."
  )
})
prediction_coverage_df <- data.frame(
  tool = names(prediction_coverage),
  n_variants = as.numeric(prediction_coverage),
  pct = round(
    as.numeric(prediction_coverage) / nrow(no_prediction_with_apogee2) * 100,
    1
  )
) |>
  dplyr::arrange(desc(n_variants))

print(prediction_coverage_df)

# Check why no predictions - variants not in apogee2 database? -----------
log_info("Checking if 'No prediction' variants are in apogee2 database...")

no_pred_variants <- unique(no_prediction_somatic$variant)
apogee2_variants <- unique(apogee2$variant)

missing_in_apogee2 <- no_pred_variants[!no_pred_variants %in% apogee2_variants]

log_info(
  "{length(missing_in_apogee2)} out of {length(no_pred_variants)} 'No prediction' variants are NOT in apogee2 database"
)

if (length(missing_in_apogee2) > 0) {
  log_info("Examples of missing variants (first 20):")
  no_prediction_somatic |>
    dplyr::filter(variant %in% missing_in_apogee2) |>
    head(20) |>
    print(width = Inf)

  log_info("Column names in no_prediction_somatic:")
  print(names(no_prediction_somatic))

  # Analyze types of variants with no prediction
  log_info("Analysis of 'No prediction' variant types:")

  no_prediction_somatic |>
    dplyr::mutate(
      is_synonymous = grepl("^\\w+:[A-Z]\\d+[A-Z]$", aachange) &
        substr(aachange, nchar(aachange), nchar(aachange)) ==
          substr(
            aachange,
            regexpr(":", aachange) + 1,
            regexpr(":", aachange) + 1
          ),
      is_stop_gain = grepl("\\*", aachange),
      is_noncoding = aachange == "non-coding"
    ) |>
    dplyr::summarise(
      n_total = n(),
      n_synonymous = sum(is_synonymous, na.rm = TRUE),
      n_stop_gain = sum(is_stop_gain, na.rm = TRUE),
      n_noncoding = sum(is_noncoding, na.rm = TRUE),
      n_other_missense = n() -
        sum(is_synonymous | is_stop_gain | is_noncoding, na.rm = TRUE)
    ) |>
    print()

  # Show examples by type
  log_info("Synonymous variants (first 10):")
  no_prediction_somatic |>
    dplyr::filter(
      grepl("^\\w+:[A-Z]\\d+[A-Z]$", aachange) &
        substr(aachange, nchar(aachange), nchar(aachange)) ==
          substr(
            aachange,
            regexpr(":", aachange) + 1,
            regexpr(":", aachange) + 1
          )
    ) |>
    dplyr::select(variant, Locus, aachange, ntchange) |>
    head(10) |>
    print()

  log_info("Non-synonymous variants:")
  no_prediction_somatic |>
    dplyr::filter(
      !grepl("^\\w+:[A-Z]\\d+[A-Z]$", aachange) |
        substr(aachange, nchar(aachange), nchar(aachange)) !=
          substr(
            aachange,
            regexpr(":", aachange) + 1,
            regexpr(":", aachange) + 1
          )
    ) |>
    dplyr::select(variant, Locus, aachange, ntchange, `Gnomad Frequency`) |>
    print(n = Inf)

  # Create pie plot of variant categories --------------------------------
  log_info("Creating pie plot of 'No prediction' variant categories...")

  variant_categories <- no_prediction_somatic |>
    dplyr::mutate(
      is_synonymous = grepl("^\\w+:[A-Z]\\d+[A-Z]$", aachange) &
        substr(aachange, nchar(aachange), nchar(aachange)) ==
          substr(
            aachange,
            regexpr(":", aachange) + 1,
            regexpr(":", aachange) + 1
          ),
      is_stop_gain = grepl("\\*", aachange),
      is_noncoding = aachange == "non-coding",
      category = case_when(
        is_synonymous ~ "Synonymous",
        is_stop_gain ~ "Stop-gain",
        is_noncoding ~ "Non-coding",
        TRUE ~ "Other"
      )
    ) |>
    dplyr::count(category) |>
    dplyr::mutate(
      percentage = round(n / sum(n) * 100, 1),
      label = glue("{category}\n(n={n}, {percentage}%)")
    )

  # Define colors for categories
  category_colors <- c(
    "Synonymous" = "#9ECAE1", # Light blue
    "Stop-gain" = "#FC9272", # Light red
    "Non-coding" = "#A1D99B" # Light green
  )

  # Create pie chart
  p_pie <- variant_categories |>
    ggplot(aes(x = "", y = n, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = category_colors) +
    geom_text(
      aes(label = glue("{n}\n({percentage}%)")),
      position = position_stack(vjust = 0.5),
      size = 5,
      fontface = "bold",
      color = "white"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30")
    ) +
    labs(
      title = "Categories of 'No Prediction' Somatic Variants",
      subtitle = glue("Total: {sum(variant_categories$n)} variants"),
      fill = "Variant Type"
    )

  # Save plot
  plot_file <- outdir / "11.02-no-prediction-variant-categories-pie.pdf"
  log_info("Saving pie plot to: {plot_file}")

  ggsave(
    filename = plot_file,
    plot = p_pie,
    width = 8,
    height = 6,
    device = "pdf"
  )

  # Also save as PNG for quick viewing
  plot_file_png <- outdir / "11.02-no-prediction-variant-categories-pie.png"
  ggsave(
    filename = plot_file_png,
    plot = p_pie,
    width = 8,
    height = 6,
    dpi = 300,
    device = "png"
  )

  log_info("Pie plot saved successfully!")
}

# For variants that ARE in apogee2, check their data
if (length(missing_in_apogee2) < length(no_pred_variants)) {
  in_apogee2 <- no_pred_variants[no_pred_variants %in% apogee2_variants]
  log_info("{length(in_apogee2)} variants are IN apogee2 database:")

  apogee2 |>
    dplyr::filter(variant %in% in_apogee2) |>
    dplyr::select(
      variant,
      Molecule_type,
      Gene_symbol,
      AA_ref,
      AA_alt,
      AA_pos,
      Functional_effect_general,
      APOGEE2,
      PolyPhen2,
      SIFT,
      AlphaMissense,
      CADD
    ) |>
    print(n = Inf)
}

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
