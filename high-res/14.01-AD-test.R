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


# Conn ---------------------------------------------------------------

conn <- db_conn(
  Sys.getenv("DUCKDB_PATH"),
  readonly = TRUE
)
tbl_ls(conn)
tbl_allvariants_cell <- dplyr::tbl(
  conn,
  "allvariants_cell"
) |>
  filter(srrid %in% ad_srrid) |>
  filter(variant_type %in% c("colorful", "black"))

# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

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
  ) -> ad_forttest

ad_forttest |>
  # head(20) |>
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
  ) -> ad_forttest_ttest


# export(
#   ad_forttest_ttest,
#   outdirnotuse / "AD" / "AD-variant-af-ttest.qs"
# )

ad_forttest_ttest <- import(
  outdirnotuse / "AD" / "AD-variant-af-ttest.qs"
)

ad_forttest_ttest |>
  select(-data) |>
  filter(variant %in% admeta_af$variant) |>
  unnest(t) |>
  dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(
    plog10p = -log10(p.value),
    est = abs(estimate),
  ) |>
  dplyr::mutate(
    rank = plog10p * est,
  ) |>
  dplyr::arrange(
    desc(rank)
  ) |>
  dplyr::rename(
    ad = "Alzheimer's Disease",
  ) |>
  dplyr::filter(
    ad >= 5,
    Healthy >= 5
  ) -> a


load_pkg(ggstatsplot)
load_pkg(ggpubr)

# ad_forttest |>
#   unnest(data) |>
#   # filter(variant %in% c("16311T>C", "263A>G", "8794C>T")) |>
#   filter(variant %in% a$variant) |>
#   ggplot(aes(x = disease, y = af)) +
#   geom_boxplot() +
#   geom_point() +
#   facet_grid(variant ~ celltype) +
#   stat_compare_means(method = "t.test")

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
