#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-25 00:47:58
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

variant_list <-
  import(
    outdirnotuse /
      "allvariants-prioritize" /
      "allvariants-prioritize-variant-annotation-with-samples-n-clean-group.qs"
  )
# Source ---------------------------------------------------------------------
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_celltype_specific_variant.R")
)
source(
  path(Sys.getenv("HIGHRESDIR"), "plot_individual_proportion.R")
)
source(
  path(Sys.getenv("HIGHRESDIR"), "00-colors.R")
)
# Conn ---------------------------------------------------------------

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

variant_annotation <- "VUS"

variant_df <- variant_list[[variant_annotation]]

variant_df |>
  left_join(
    ALLVARIANTS |>
      left_join(
        METAFULL |>
          select(srrid, disease),
        by = c("srrid")
      ) |>
      mutate(
        variant_type = factor(variant_type, c("hete", "homo")),
        disease = factor(
          disease,
          levels = names(color_disease)
        )
      ) |>
      arrange(variant_type, disease) |>
      nest(.by = "variant", .key = "meta"),
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
      prediction_class = prediction_class_new,
      Disease = Disease,
      mc.cores = 10,
      SIMPLIFY = FALSE
    )
  ) -> admeta_af_plot


admeta_af_plot |>
  filter(n_hete > 0) |>
  # head(2) |>
  mutate(
    saveplot = mapply(
      FUN = \(.variant, .plot, .prediction_class, .disease) {
        tryCatch(
          expr = {
            safe_variant <- sanitize_path_component(.variant)
            # safe_adspecific <- sanitize_path_component(.adspecific)
            safe_prediction_class <- sanitize_path_component(.prediction_class)
            safe_disease <- sanitize_path_component(.disease)

            saveplot(
              plot = .plot,
              filename = outdirnotuse /
                "allvariants-prioritize" /
                "candidates" /
                # safe_adspecific /
                safe_prediction_class /
                safe_disease /
                glue(
                  "{safe_variant}-{safe_prediction_class}-{safe_disease}-CELLTYPE-SPECIFIC-HIST-PLOT.pdf"
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
      .prediction_class = prediction_class_new,
      .disease = Disease
    )
  )

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
