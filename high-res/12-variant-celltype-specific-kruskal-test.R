#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-31 15:42:52
# @DESCRIPTION: this script is used for ...

# Library -----------------------------------------------------------------

load_pkg(jutils)

# args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
GetoptLong.options(help_style = "two-column")
VERSION = "v0.0.1"

# default: default value specified here.

verbose = TRUE

GetoptLong("verbose!", "print messages")


logger::log_threshold(logger::TRACE)
logger::log_layout(logger::layout_glue_colors)

# header ------------------------------------------------------------------

# load data ---------------------------------------------------------------
outdir <- path("/home/liuc9/github/scMOCHA-data/analysis/zzz/MANUSCRIPTFIGURES")
cleandatadir <- path("/home/liuc9/github/scMOCHA-data/data/zzz/clean-data")

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
)

HOMO_HETE_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type %in% c("homo", "hete"))
dotenv(".env")
# load conn ---------------------------------------------------------------

conn <- conn_db(
  Sys.getenv("DUCKDB_PATH"),
  read_only = TRUE
)
DBI::dbListTables(conn)
tbl_allvariants_cell <- dplyr::tbl(
  conn,
  "allvariants_cell"
)


# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------
fn_kruskal_test <- function(.gseid, .srrid, .variant) {
  # .gseid <- "GSE226602"
  # .srrid <- "GSM7080053"
  # .variant <- "11251A>G"
  log_info(glue::glue("Processing {.gseid} - {.srrid} - {.variant}"))

  tbl_allvariants_cell |>
    dplyr::filter(
      srrid == .srrid,
      variant == .variant,
      variant_type %in% c("colorful", "black")
    ) |>
    as.data.table() |>
    # dplyr::mutate(
    #   af = ifelse(
    #     af < 0.01,
    #     NA_real_,
    #     af
    #   )
    # ) |>
    dplyr::filter(celltype != "other") -> dt_variant_celltype_af

  tryCatch(
    {
      kruskal.test(
        af ~ celltype,
        data = dt_variant_celltype_af
      ) |>
        broom::tidy()
    },
    error = \(e) {
      return(tibble::tibble(
        statistic = NA_real_,
        p.value = NA_real_,
        method = "Kruskal-Wallis test",
        parameter = NA_real_
      ))
    }
  )
}
# body --------------------------------------------------------------------

HOMO_HETE_VARIANTS |>
  # head(200) |>
  dplyr::mutate(
    kruskal_test = parallel::mcmapply(
      FUN = fn_kruskal_test,
      .gseid = gseid,
      .srrid = srrid,
      .variant = variant,
      SIMPLIFY = FALSE,
      mc.cores = 20
    )
  ) |>
  tidyr::unnest(kruskal_test) -> HOMO_HETE_VARIANTS_KRUSKAL


{
  export(
    HOMO_HETE_VARIANTS_KRUSKAL,
    outdir /
      "VARIANT-KRUSKAL-WALLIS-TEST.xlsx"
  )
}
close_all_db()
# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
