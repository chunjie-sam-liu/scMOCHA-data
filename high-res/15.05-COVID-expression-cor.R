#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-03-23 00:00:00
# @DESCRIPTION: COVID-19 top-variant expression correlation analysis

# Reproducibility ----------------------------------------------------------
set.seed(1)

# Library -----------------------------------------------------------------

suppressMessages({
  load_pkg(jutils)
})

# Args --------------------------------------------------------------------

VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

nthread = 30
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

# Load data ----------------------------------------------------------------

dotenv()
suppressMessages({
  conflicted::conflict_prefer("filter", "dplyr")
})

outdir <- path(Sys.getenv("OUTDIR"))
outdirnotuse <- path(Sys.getenv("OUTDIRNOTUSE"))
covid_outdir <- outdirnotuse / "COVID19"
corr_outdir <- covid_outdir / "corr"
expr_path <- Sys.getenv(
  "SCMOCHA_EXPR_PATH",
  unset = "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/EXPR/gse_srrid_celltype_gene_expr.fst"
)
dir_create(corr_outdir)

allvariants <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  filter(variant_type %in% c("hete", "homo"))
metafull <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")

covid_meta <- metafull |>
  filter(disease %in% c("Healthy", "COVID-19")) |>
  select(gseid, srrid, disease) |>
  mutate(disease = factor(disease, levels = c("Healthy", "COVID-19")))

cluster_variant <- import(outdir / "ALLVARIANT-ALLSAMPLES-CLUSTERAF.qs")

expr <- import(expr_path) |>
  inner_join(covid_meta |> select(srrid, disease), by = "srrid") |>
  pivot_longer(
    cols = -c(gseid, srrid, genename, disease),
    names_to = "celltype",
    values_to = "expr"
  ) |>
  nest(.by = c(genename, celltype), .key = "expr")

variants <- import(covid_outdir / "COVID19-variant-top-ttest-cluster-variants.qs2")

# Functions ----------------------------------------------------------------

normalize_qs2_output <- function(.path_qs2) {
  .path_qs2 <- as.character(.path_qs2)
  .path_qs <- sub("\\.qs2$", ".qs", .path_qs2)
  if (file.exists(.path_qs)) {
    file.copy(.path_qs, .path_qs2, overwrite = TRUE)
    invisible(file.remove(.path_qs))
  }
}

fn_cor_test_variant <- function(.variant) {
  csv_path <- corr_outdir / glue("covid-celltype-variant-af-{.variant}-corr.csv")
  qs2_path <- corr_outdir / glue("covid-celltype-variant-af-{.variant}-corr.qs2")

  if (file_exists(csv_path) && file_exists(qs2_path)) {
    log_info("Skipping {.variant}; correlation outputs already exist")
    return(invisible(NULL))
  }

  log_info("Processing COVID correlation for {.variant}")

  af_nested <- cluster_variant |>
    select(celltype, gseid, srrid, af = all_of(.variant)) |>
    filter(srrid %in% covid_meta$srrid) |>
    nest(.by = celltype, .key = "af")

  expr |>
    left_join(af_nested, by = "celltype") -> expr_variant

  expr_variant |>
    mutate(
      corr_results = pbmcmapply(
        function(.expr, .af) {
          tryCatch(
            {
              .expr |>
                left_join(.af, by = c("gseid", "srrid")) -> expr_af

              expr_af |>
                filter(!is.na(af), !is.na(expr)) -> expr_af_cc

              expr_af_cc |>
                count(disease, .drop = FALSE) |>
                pivot_wider(
                  names_from = disease,
                  values_from = n,
                  values_fill = 0
                ) -> n_disease

              if (n_disease$Healthy < 10 || n_disease$`COVID-19` < 10) {
                return(NULL)
              }

              cor.test(~ af + expr, data = expr_af_cc, method = "pearson") |>
                broom::tidy() |>
                select(corr = estimate, pval = p.value) |>
                bind_cols(n_disease)
            },
            error = function(e) {
              log_warn("Correlation failed for {.variant}: {conditionMessage(e)}")
              NULL
            }
          )
        },
        .expr = expr,
        .af = af,
        mc.cores = nthread,
        mc.preschedule = TRUE,
        SIMPLIFY = FALSE
      )
    ) -> tmp

  tmp |>
    mutate(
      corr_results = pbmclapply(
        corr_results,
        function(.x) {
          if (is.null(.x) || inherits(.x, "try-error")) NULL else .x
        },
        mc.cores = nthread,
        mc.preschedule = TRUE
      )
    ) |>
    filter(!purrr::map_lgl(corr_results, is.null)) |>
    unnest(cols = corr_results) |>
    group_by(celltype) |>
    mutate(adj_pval = p.adjust(pval, method = "BH")) |>
    ungroup() -> expr_variant_corr

  export(expr_variant_corr |> select(-c(expr, af)), csv_path)
  export(expr_variant_corr, qs2_path)
  normalize_qs2_output(qs2_path)
}

# Main ---------------------------------------------------------------------

variants |>
  walk(fn_cor_test_variant)

if (isTRUE(verbose)) {
  sessionInfo()
}
