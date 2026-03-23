#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: Cluster-level COVID-19 vs Healthy variant AF t-tests

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library ------------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args ---------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)

# Logger -------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}

# Load data ----------------------------------------------------------------

load_pkg(jutils)
dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"

dir.create(covid_outdir, recursive = TRUE, showWarnings = FALSE)

cluster_variant <- import(
  outdir / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
)
allvariants <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
metafull <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
variant_nsamples_wide_annotated <- import(
  covid_outdir / "COVID19-variant-sample-counts-annotated.xlsx"
)

# Main ---------------------------------------------------------------------

metafull |>
  dplyr::filter(disease %in% c("Healthy", "COVID-19")) |>
  select(gseid, srrid, Chemistry, disease) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "COVID-19"))
  ) -> covid_meta

covid_meta |>
  left_join(
    allvariants,
    by = c("gseid", "srrid")
  ) |>
  select(
    -c(Chemistry, Haplogroup, Verbose_haplogroup)
  ) |>
  mutate(
    disease = factor(disease, levels = c("Healthy", "COVID-19"))
  ) -> covid_meta_af

variant_nsamples_wide_annotated |>
  filter(
    Healthy >= 5,
    `COVID-19` >= 5
  ) |>
  pull(variant) |>
  unique() -> covid_variants_test

covid_srrid <- covid_meta_af$srrid |> unique()
covid_variant <- intersect(covid_meta_af$variant |> unique(), covid_variants_test)

log_info("Testing {length(covid_variant)} variants across cluster-level cell types")

cluster_variant |>
  filter(srrid %in% covid_srrid) |>
  pivot_longer(
    cols = -c(gseid, srrid, celltype),
    names_to = "variant",
    values_to = "af"
  ) |>
  filter(variant %in% covid_variant) |>
  inner_join(
    covid_meta |>
      select(gseid, srrid, disease),
    by = c("gseid", "srrid")
  ) |>
  nest(
    .by = c(variant, celltype)
  ) -> covid_forttest_cluster

covid_forttest_cluster |>
  mutate(
    t = pbmclapply(
      X = data,
      FUN = function(.x) {
        .x |>
          dplyr::count(disease, .drop = FALSE) |>
          tidyr::pivot_wider(
            names_from = disease,
            values_from = n,
            values_fill = 0
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
              dplyr::bind_cols(.xx)
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
                `COVID-19` = NA_integer_
              )
            )
          }
        )
      },
      mc.cores = nthread
    )
  ) -> covid_forttest_cluster_ttest

export(
  covid_forttest_cluster_ttest,
  covid_outdir / "COVID19-variant-af-ttest-cluster.qs2"
)

output_qs2 <- as.character(covid_outdir / "COVID19-variant-af-ttest-cluster.qs2")
output_qs <- as.character(covid_outdir / "COVID19-variant-af-ttest-cluster.qs")

if (file.exists(output_qs)) {
  file.copy(output_qs, output_qs2, overwrite = TRUE)
  invisible(file.remove(output_qs))
}

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
