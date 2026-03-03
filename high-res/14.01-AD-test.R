#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-02-23 11:27:36
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

source("high-res/00-colors.R")

# Load data ---------------------------------------------------------------
load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
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

ad_srrid <- admeta_af$srrid |> unique()
ad_variant <- admeta_af$variant |> unique()

# Conn ---------------------------------------------------------------

conn <- db_conn(
  Sys.getenv("DUCKDB_PATH"),
  readonly = TRUE
)
tbl_ls(conn)
dt_allvariants_cell <- dplyr::tbl(
  conn,
  "allvariants_cell"
) |>
  filter(srrid %in% ad_srrid) |>
  filter(variant %in% ad_variant) |>
  filter(variant_type %in% c("colorful", "black")) |>
  as.data.table() |>
  inner_join(
    admeta |>
      select(gseid, srrid, disease),
    by = c("gseid", "srrid")
  )

db_disconn()
# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

dt_allvariants_cell |>
  select(gseid, srrid, celltype, variant, af, disease, barcode) |>
  nest(
    .by = c(variant, celltype)
  ) -> ad_forttest_cell

cluster_variant |>
  filter(srrid %in% ad_srrid) |>
  pivot_longer(
    cols = -c(gseid, srrid, celltype),
    names_to = "variant",
    values_to = "af"
  ) |>
  inner_join(
    admeta |>
      select(gseid, srrid, disease),
    by = c("gseid", "srrid")
  ) |>
  nest(
    .by = c(variant, celltype)
  ) -> ad_forttest_cluster

ad_forttest_cluster |>
  mutate(
    t = pbmclapply(
      X = data,
      FUN = function(.x) {
        .x |>
          dplyr::count(disease) |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n
          ) -> .xx

        tryCatch(
          expr = {
            t.test(
              af ~ disease,
              data = .x,
              var.equal = TRUE
            ) |>
              broom::tidy() |>
              dplyr::select(
                estimate,
                estimate1,
                estimate2,
                p.value,
                conf.low,
                conf.high
              ) |>
              dplyr::bind_cols(
                .xx
              )
          },
          error = function(e) {
            message("Error: ", conditionMessage(e))
            return(
              tibble::tibble(
                estimate = NA_real_,
                estimate1 = NA_real_,
                estimate2 = NA_real_,
                p.value = NA_real_,
                conf.low = NA_real_,
                conf.high = NA_real_,
                Healthy = NA_integer_,
                `Alzheimer's Disease` = NA_integer_
              )
            )
          }
        )
      },
      mc.cores = 20
    )
  ) -> ad_forttest_cluster_ttest

export(
  ad_forttest_cluster_ttest,
  outdirnotuse / "AD" / "AD-variant-af-ttest-cluster.qs"
)

ad_forttest_cell |>
  mutate(
    t = pbmclapply(
      X = data,
      FUN = function(.x) {
        .x |>
          dplyr::count(disease) |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n
          ) -> .xx

        tryCatch(
          expr = {
            t.test(
              af ~ disease,
              data = .x,
              var.equal = TRUE
            ) |>
              broom::tidy() |>
              dplyr::select(
                estimate,
                estimate1,
                estimate2,
                p.value,
                conf.low,
                conf.high
              ) |>
              dplyr::bind_cols(
                .xx
              )
          },
          error = function(e) {
            message("Error: ", conditionMessage(e))
            return(
              tibble::tibble(
                estimate = NA_real_,
                estimate1 = NA_real_,
                estimate2 = NA_real_,
                p.value = NA_real_,
                conf.low = NA_real_,
                conf.high = NA_real_,
                Healthy = NA_integer_,
                `Alzheimer's Disease` = NA_integer_
              )
            )
          }
        )
      },
      mc.cores = 20
    )
  ) -> ad_forttest_cell_ttest


export(
  ad_forttest_cell_ttest,
  outdirnotuse / "AD" / "AD-variant-af-ttest-cell.qs"
)


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
