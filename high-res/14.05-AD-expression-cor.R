#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-03 00:59:53
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
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
METAFULL |>
  dplyr::filter(
    disease %in% c("Healthy", "Alzheimer's Disease"),
    Chemistry == "SC5P-PE"
  ) |>
  dplyr::select(gseid, srrid, Chemistry, disease) -> admeta

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

cluster_variant <- import(
  path(Sys.getenv("OUTDIR")) / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs"
)

expr <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
) |>
  dplyr::inner_join(
    admeta |> select(srrid, disease),
    by = c("srrid")
  ) |>
  tidyr::pivot_longer(
    cols = -c(gseid, srrid, genename, disease),
    names_to = "celltype",
    values_to = "expr"
  ) |>
  tidyr::nest(
    .by = c(genename, celltype),
    .key = "expr"
  )

variants <- c(
  "13592C>T",
  "5031G>T",
  "8362T>G"
)

variants <- import(
  outdirnotuse / "AD" / "AD-variant-top-ttest-cluster-variants.qs"
)

# Conn ---------------------------------------------------------------

# Function ----------------------------------------------------------------

fn_cor_test_variant <- function(.variant) {
  .v <- cluster_variant |>
    select(celltype, gseid, srrid, af = all_of(.variant)) |>
    filter(srrid %in% admeta$srrid) |>
    nest(.by = celltype, .key = "af")

  expr |>
    dplyr::left_join(
      .v,
      by = c("celltype")
    ) -> .expr_v

  .expr_v |>
    # head(20) |>
    dplyr::mutate(
      corr_results = pbmcmapply(
        function(.expr, .af) {
          # .expr <- .expr_v$expr[[1]]
          # .af <- .expr_v$af[[1]]
          tryCatch(
            {
              .expr |>
                dplyr::left_join(
                  .af,
                  by = c("gseid", "srrid")
                ) -> .expr_af

              .expr_af |>
                dplyr::count(disease) |>
                tidyr::pivot_wider(
                  names_from = disease,
                  values_from = n
                ) -> .n_disease

              if (sum(.n_disease < 10) > 0) {
                return(NULL)
              }

              cor.test(
                ~ af + expr,
                data = .expr_af,
                method = "pearson"
              ) -> .corr

              .corr |>
                broom::tidy() |>
                dplyr::select(
                  corr = estimate,
                  pval = p.value,
                ) |>
                dplyr::bind_cols(
                  .n_disease
                )
            },
            error = function(e) {
              return(NULL)
            }
          )
        },
        .expr = expr,
        .af = af,
        mc.cores = 30,
        mc.preschedule = TRUE,
        SIMPLIFY = FALSE
      )
    ) -> .tmp

  .tmp |>
    dplyr::mutate(
      corr_results = pbmclapply(
        corr_results,
        function(.x) {
          if (is.null(.x) || inherits(.x, "try-error")) NULL else .x
        },
        mc.cores = 30,
        mc.preschedule = TRUE
      )
    ) |>
    dplyr::filter(!purrr::map_lgl(corr_results, is.null)) |>
    tidyr::unnest(cols = corr_results) -> .expr_v_corr

  export(
    .expr_v_corr |> select(-c(expr, af)),
    outdirnotuse /
      "AD" /
      "corr" /
      "ad-celltype-variant-af-{.variant}-corr.csv" |> glue()
  )
  export(
    .expr_v_corr,
    outdirnotuse /
      "AD" /
      "corr" /
      "ad-celltype-variant-af-{.variant}-corr.qs" |> glue()
  )
}


# Main --------------------------------------------------------------------

variants[-1] |>
  purrr::walk(
    fn_cor_test_variant
  )


# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
