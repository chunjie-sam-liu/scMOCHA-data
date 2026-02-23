#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-23 13:01:02
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
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

variant_annotation <- import(
  path(Sys.getenv("OUTDIR")) /
    "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.xlsx"
)
#

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(
  Sys.getenv("OUTDIRNOTUSE")
)
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")

METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease")
  ) |>
  filter(Chemistry == "SC5P-PE") |>
  select(gseid, srrid, Chemistry, disease) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta

ad_srrid <- admeta$srrid

admeta |>
  left_join(
    ALLVARIANTS,
    by = c("gseid", "srrid")
  ) |>
  select(
    -c(Chemistry, Haplogroup, Verbose_haplogroup)
  ) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "Alzheimer's Disease"))
  ) -> admeta_af

ad_forttest_ttest <- import(
  outdirnotuse / "AD" / "AD-variant-af-ttest.qs"
)

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

variant_annotation |>
  filter(
    prediction_class %in% c("pathogenic", "likely_pathogenic") | !is.na(Disease)
  ) |>
  filter(variant %in% admeta_af$variant)


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
