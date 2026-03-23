#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 01:51:14
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
  outdirnotuse / "COVID19" / "COVID19-variant-sample-counts-annotated.xlsx"
)

# Source ---------------------------------------------------------------------
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_celltype_specific_variant.R")
)
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_individual_proportion.R")
)
# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

METAFULL |>
  filter(
    disease %in% c("Healthy", "COVID-19"),
    # Chemistry == "SC5P-PE"
  ) |>
  select(gseid, srrid, Chemistry, disease) -> disease_meta

disease_meta |>
  dplyr::left_join(ALLVARIANTS, by = c("gseid", "srrid")) |>
  dplyr::select(-c(Chemistry, Haplogroup, Verbose_haplogroup)) |>
  dplyr::mutate(
    disease = factor(disease, levels = c("Healthy", "COVID-19"))
  ) -> meta_af


variant_nsamples_wide_annotated |>
  filter(
    `Healthy` == 0
  ) |>
  filter(
    `COVID-19` > 1
  ) |>
  pull(variant) -> variants_disease

variant_nsamples_wide_annotated |>
  filter(
    `COVID-19` > 5,
    `Healthy` > 5
  ) |>
  pull(variant) -> variants_test

thevariants <- c(
  variants_disease,
  variants_test
)

export(thevariants, outdirnotuse / "COVID19" / "candidates" / "thevariants.qs")


meta_af |>
  filter(variant %in% thevariants) |>
  select(gseid, srrid, disease, variant_type, variant) |>
  nest(.by = variant, .key = "meta") |>
  mutate(
    adspecific = if_else(
      variant %in% variants_disease,
      "COVID-19-specific",
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
    plot = pbmcmapply(
      FUN = \(.x, .y, aachange, prediction_class, Disease) {
        # .x <- forplots$variant[1]
        # .y <- forplots$meta[[1]]
        # fn_plot_hist(thevariant, thegseid, thesrrid, subtitle = "cj")
        .y |>
          arrange(desc(disease), variant_type) |>
          mutate(
            info = glue(
              "{disease} - {variant_type} - {aachange} - {prediction_class} - {Disease}"
            )
          ) |>
          # head(2) |>
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
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) -> admeta_af_plot


export(
  admeta_af_plot,
  outdirnotuse / "COVID19" / "candidates" / "admeta_af_plot.qs"
)
admeta_af_plot |>
  head(2) |>
  mutate(
    saveplot = mapply(
      FUN = \(.variant, .plot, .adspecific, .prediction_class, .disease) {
        tryCatch(
          expr = {
            saveplot(
              plot = .plot,
              file = outdirnotuse /
                "COVID19" /
                "candidates" /
                .adspecific /
                glue(
                  "COVID19-{.variant}-{.adspecific}-{.prediction_class}-{.disease}-CELLTYPE-SPECIFIC-HIST-PLOT.pdf"
                ),
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


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
