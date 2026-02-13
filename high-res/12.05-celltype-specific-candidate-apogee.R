#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-13 05:43:26
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


outdir <- path(Sys.getenv("OUTDIR"))

color_celltype <- c(
  "homo hetero" = "blue",
  "Bimodal" = "red",
  "Mono" = "#A6D854FF",
  "B" = "#66C2A5FF",
  "T cell" = "#E5C494FF",
  "CD8 T" = "#8DA0CBFF",
  "NK" = "#FFD92FFF",
  "other" = "black"
)

variant_annotation <- import(
  outdir / "VARIANT-ANNOTATION-TABLE-APOGEE2.xlsx"
)

thevariants <- tibble(
  variant = c(
    "709G>A", # homo hetero
    "14905G>A", # Bimodal
    "8138A>G", # Bimodal
    "2011G>A", # Bimodal
    "7751T>C", # Bimodal
    "7609T>C", # Bimodal
    "4813T>C", # Mono
    "7159T>C", # B cell
    "7833T>C", # T cell
    "10500G>A", # T cell
    "10097A>G", # T cell
    "8005T>C", # CD8 T
    "7850G>A", # CD8 T
    "9033A>G", # NK
    "7757G>A", # NK
    "9390A>G", # NK
    "6374T>C", # NK
    "10236A>G", # NK
    "1474G>A", # NK
    "9609T>C", # NK
    "2636G>A", # NK
    "15612G>A", # NK
    "2343G>A", # NK
    "7837T>C", # NK
    "6928T>C", # NK
    "2666T>C" # NK
  ),
  type = c(
    "homo hetero", # homo hetero
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal",
    "Bimodal", # Bimodal
    "Mono", # Mono
    "B", # B cell
    "T cell",
    "T cell",
    "T cell", # T cell
    "CD8 T",
    "CD8 T", # CD8 T
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK",
    "NK" # NK
  )
)

ALLVARIANTS_TEST <- import(
  outdir / "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
)

ALLVARIANTS_TEST_SIG <- ALLVARIANTS_TEST |>
  dplyr::filter(p.value < 0.05, statistic > 20)


thevariants |>
  left_join(variant_annotation, by = c("variant" = "Variant")) |>
  mutate(Locus = gsub(",Humanin", "", Locus)) -> thevariants_annotation


thevariants_annotation |>
  filter(!is.na(`Prediction Class`)) |>
  count(`Prediction Class`)


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

# variant_annotation <- import(
#   outdir / "VARIANT-ANNOTATION-TABLE.xlsx"
# )

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
thevariants_annotation |> glimpse()
thevariants_annotation |> filter(!is.na(Apogee2))

thevariants_annotation |>
  filter(
    !is.na(Disease) |
      !`Prediction Class` %in% c("Benign", "Likely benign", "No prediction")
  ) |>
  filter(!`Aachange New` %in% c("rRNA", "tRNA")) |>
  glimpse()


thevariants |>
  left_join(variant_annotation, by = c("variant" = "Variant")) |>
  mutate(Locus = gsub(",Humanin", "", Locus)) -> thevariants_annotation


thevariants_annotation |>
  filter(!is.na(`Prediction Class`)) |>
  count(`Prediction Class`)

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------
thevariants_annotation |> glimpse()
thevariants_annotation |> filter(!is.na(Apogee2))


ALLVARIANTS_TEST |>
  dplyr::mutate(
    log10p = -log10(p.value),
    log10p = ifelse(is.infinite(log10p), 300, log10p)
  ) |>
  select(gseid, srrid, variant, statistic, log10p) |>
  dplyr::filter(!is.na(log10p)) |>
  left_join(variant_annotation, by = c("variant" = "Variant")) |>
  filter(
    !is.na(Disease) |
      grepl(
        "pathogenic",
        `Prediction Class`,
        ignore.case = TRUE
      )
  ) |>
  filter(!`Aachange Group` %in% c("rRNA", "tRNA")) |>
  arrange(-log10p) -> candidate_variants_apogee2

candidate_variants_apogee2 |>
  filter(
    variant %in% SOMATIC_VARIANTS$variant
  ) |>
  export(
    "/home/liuc9/github/scMOCHA-data/high-res-MANUSCRIPTFIGURES-notuse/12.celltype-specific-candidate-apogee2.xlsx"
  )


candidate_variants_apogee2 |>
  filter(
    variant %in% SOMATIC_VARIANTS$variant
  ) |>
  glimpse()

# Main --------------------------------------------------------------------

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
