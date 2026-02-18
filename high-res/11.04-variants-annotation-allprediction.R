#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-17 14:39:48
# @DESCRIPTION: Explore and annotate somatic variants with multiple prediction tools
# @VERSION: v0.1.0
#
# Prediction Software Coverage:
# - APOGEE2 (MitImpact): Mitochondrial-specific pathogenicity predictor (32 scores)
# - nAPOGEE: Non-coding variant pathogenicity predictor
# - AlphaMissense: DeepMind AI-based missense variant prediction
# - CADD: Combined Annotation Dependent Depletion
# - PolyPhen2: Polymorphism Phenotyping v2
# - SIFT/SIFT4G: Sort Intolerant From Tolerant
# - MutationTaster: Disease-causing potential
# - FATHMM: Functional Analysis through HMM
# - Additional tools: PROVEAN, Condel, CAROL, MitoTip, etc.
#
# Key Categories:
# 1. Evolutionary Conservation: SIFT, PROVEAN
# 2. Structural-based: PolyPhen2, MutationAssessor
# 3. ML Ensemble: APOGEE2, nAPOGEE, AlphaMissense, Meta-SNP
# 4. Functional/Phenotype: CADD, VEST, MutationTaster, FATHMM
# 5. Mitochondria-specific: MitoTip, SNPDryad, Mitoclass1

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
# source("high-res/00-colors.R")  # Not needed for variant annotation

# Load data ---------------------------------------------------------------
load_pkg(jutils)
# dotenv(".env")  # .env is in parent directory
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

# Set paths
basedir <- path("/liulab/chunjie/data/scMOCHA")
outdir <- path(basedir, "analysis/zzz/MANUSCRIPTFIGURES")
outdirnotuse <- path(basedir, "analysis/zzz/MANUSCRIPTFIGURES-notuse")

# Check if running from high-res directory
if (basename(getwd()) == "high-res") {
  outdir <- path("../high-res-MANUSCRIPTFIGURES")
  outdirnotuse <- path("../high-res-MANUSCRIPTFIGURES-notuse")
}


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

# Merge somatic variants with predictions
somatic_with_pred <- SOMATIC_VARIANTS |>
  left_join(
    apogee2 |>
      select(
        variant,
        Gene_symbol,
        Functional_effect_general,
        PolyPhen2,
        PolyPhen2_score,
        SIFT,
        SIFT_score,
        APOGEE2,
        APOGEE2_score,
        APOGEE2_probability,
        AlphaMissense,
        AlphaMissense_score,
        CADD_phred_score,
        MutationTaster,
        FATHMM,
        FATHMM_score,
        MitoTip_count,
        Condel,
        Condel_score
      ),
    by = "variant",
    relationship = "many-to-many"
  ) |>
  left_join(
    apogee2_n |>
      select(
        variant,
        nAPOGEE_score,
        nAPOGEE_posterior_probability,
        pathogenicity_assessment
      ),
    by = "variant",
    suffix = c("_mitimpact", "_napogee"),
    relationship = "many-to-many"
  )

# Summary statistics
log_info("Total somatic variants: {nrow(SOMATIC_VARIANTS)}")
log_info(
  "With APOGEE2 predictions: {sum(!is.na(somatic_with_pred$APOGEE2_score))}"
)
log_info(
  "With nAPOGEE predictions: {sum(!is.na(somatic_with_pred$nAPOGEE_score))}"
)
log_info(
  "With AlphaMissense: {sum(!is.na(somatic_with_pred$AlphaMissense_score))}"
)
log_info("With CADD: {sum(!is.na(somatic_with_pred$CADD_phred_score))}")

# APOGEE2 classification distribution
apogee2_dist <- somatic_with_pred |>
  filter(!is.na(APOGEE2)) |>
  count(APOGEE2) |>
  arrange(desc(n))

log_info("APOGEE2 classification distribution:")
print(apogee2_dist)

# Top pathogenic variants
top_pathogenic <- somatic_with_pred |>
  filter(APOGEE2 %in% c("Pathogenic", "Likely-pathogenic")) |>
  arrange(desc(APOGEE2_score)) |>
  select(
    variant,
    Gene_symbol,
    APOGEE2,
    APOGEE2_score,
    APOGEE2_probability,
    AlphaMissense,
    AlphaMissense_score,
    CADD_phred_score,
    PolyPhen2,
    SIFT
  ) |>
  distinct(variant, .keep_all = TRUE)

log_info("Found {nrow(top_pathogenic)} pathogenic/likely-pathogenic variants")

# Gene distribution
gene_dist <- top_pathogenic |>
  count(Gene_symbol, APOGEE2) |>
  arrange(desc(n))

log_info("Genes with pathogenic variants:")
print(gene_dist)

# Save  --------------------------------------------------------------
export(somatic_with_pred, outdir / "SOMATIC-VARIANTS-WITH-PREDICTIONS.xlsx")
log_success(
  "Saved merged predictions to: SOMATIC-VARIANTS-WITH-PREDICTIONS.xlsx"
)

export(top_pathogenic, outdir / "SOMATIC-VARIANTS-PATHOGENIC-TOP.xlsx")
log_success(
  "Saved top pathogenic variants to: SOMATIC-VARIANTS-PATHOGENIC-TOP.xlsx"
)

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
