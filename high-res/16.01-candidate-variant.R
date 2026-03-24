#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 10:20:00
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

variant_nsamples_wide_annotated <- import(
  outdirnotuse /
    "Other-disease" /
    "Other-disease-variant-sample-counts-annotated.xlsx"
)

# Source ------------------------------------------------------------------
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_celltype_specific_variant.R")
)
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_individual_proportion.R")
)

# Function ----------------------------------------------------------------
sanitize_path_component <- function(x, missing = "NA") {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || trimws(x) == "") {
    return(missing)
  }

  x <- gsub("[/\\\\]+", "-", x)
  x <- gsub("[:*?\"<>|]", "-", x)
  trimws(x)
}

# Main --------------------------------------------------------------------
METAFULL |>
  mutate(
    disease = dplyr::case_when(
      disease == "Healthy" ~ "Healthy",
      !is.na(disease) & !disease %in% c("Alzheimer's Disease", "COVID-19", "Unknown") ~
        "Other disease",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(disease)) |>
  select(gseid, srrid, Chemistry, disease) -> disease_meta

disease_meta |>
  dplyr::left_join(ALLVARIANTS, by = c("gseid", "srrid")) |>
  dplyr::select(-c(Chemistry, Haplogroup, Verbose_haplogroup)) |>
  dplyr::mutate(
    disease = factor(disease, levels = c("Healthy", "Other disease"))
  ) -> meta_af

variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  ) |>
  filter(
    `Other disease` > 1
  ) |>
  pull(variant) -> variants_disease

variant_nsamples_wide_annotated |>
  filter(
    `Other disease` > 5,
    `Healthy` > 5
  ) |>
  pull(variant) -> variants_test

thevariants <- c(
  variants_disease,
  variants_test
)

export(
  thevariants,
  outdirnotuse / "Other-disease" / "candidates" / "thevariants.qs"
)

meta_af |>
  filter(variant %in% thevariants) |>
  select(gseid, srrid, disease, variant_type, variant) |>
  nest(.by = variant, .key = "meta") |>
  mutate(
    adspecific = if_else(
      variant %in% variants_disease,
      "Other-disease-specific",
      "Test"
    ),
  ) |>
  left_join(
    variant_nsamples_wide_annotated |>
      select(variant, aachange, prediction_class, Disease),
    by = "variant"
  ) -> forplots

forplots |>
  mutate(
    plot = mapply(
      FUN = \(.x, .y, aachange, prediction_class, Disease) {
        .y |>
          arrange(desc(disease), variant_type) |>
          mutate(
            info = glue(
              "{disease} - {variant_type} - {aachange} - {prediction_class} - {Disease}"
            )
          ) |>
          mutate(
            plot = pmap(
              list(
                .thegseid = gseid,
                .thesrrid = srrid,
                .subtitle = info
              ),
              .f = \(.thegseid, .thesrrid, .subtitle) {
                fn_plot_hist(
                  thevariant = .x,
                  thegseid = .thegseid,
                  thesrrid = .thesrrid,
                  subtitle = .subtitle
                )
              }
            )
          ) -> plots
        plots$plot
      },
      .x = variant,
      .y = meta,
      aachange = aachange,
      prediction_class = prediction_class,
      Disease = Disease,
      SIMPLIFY = FALSE
    )
  ) -> othermeta_af_plot

lobstr::obj_size(othermeta_af_plot)
export(
  othermeta_af_plot,
  outdirnotuse / "Other-disease" / "candidates" / "othermeta_af_plot.rds"
)

othermeta_af_plot |>
  mutate(
    saveplot = mapply(
      FUN = \(.variant, .plot, .adspecific, .prediction_class, .disease) {
        tryCatch(
          expr = {
            safe_variant <- sanitize_path_component(.variant)
            safe_adspecific <- sanitize_path_component(.adspecific)
            safe_prediction_class <- sanitize_path_component(.prediction_class)
            safe_disease <- sanitize_path_component(.disease)

            saveplot(
              plot = .plot,
              filename = outdirnotuse /
                "Other-disease" /
                "candidates" /
                safe_adspecific /
                safe_prediction_class /
                safe_disease /
                glue(
                  "Other-disease-{safe_variant}-{safe_adspecific}-{safe_prediction_class}-{safe_disease}-CELLTYPE-SPECIFIC-HIST-PLOT.pdf"
                ),
              device = "pdf",
              width = 20,
              height = 10,
              create.dir = TRUE
            )
          },
          error = function(e) {
            log_error(
              glue(
                "Error saving plot for variant {.variant}: {e$message}"
              )
            )
          }
        )
      },
      .variant = variant,
      .plot = plot,
      .adspecific = adspecific,
      .prediction_class = prediction_class,
      .disease = Disease
    )
  )

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
