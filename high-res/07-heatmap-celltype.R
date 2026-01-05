#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-12-16 12:22:30
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

allvariants <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/high-res/MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
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


allvariants |>
  dplyr::filter(variant_type %in% c("homo", "hete")) -> allvariants_homohete

allpos <- unique(allvariants_homohete$start) |>
  sort() |>
  as.character()

thevariants <- unique(allvariants_homohete$variant) |>
  sort() |>
  as.character()
allsrrid <- unique(allvariants_homohete$srrid) |> as.character()
dotenv(".env")
# load conn ---------------------------------------------------------------
dotenv(".env")
conn <- DBI::dbConnect(
  duckdb::duckdb(),
  Sys.getenv("DUCKDB_PATH_COV"),
  read_only = TRUE
)

conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  Sys.getenv("DUCKDB_PATH"),
  read_only = TRUE
)


DBI::dbListTables(conn)
DBI::dbListTables(conn_all_hetero_af)
tbl_covall <- dplyr::tbl(conn, "covall")
tbl_barcode <- dplyr::tbl(conn_all_hetero_af, "barcode")
colnames(tbl_covall) |> head()

tbl_covall |>
  dplyr::filter(srrid %in% allsrrid) -> tbl_covall_allsrrid

tbl_barcode |>
  dplyr::filter(srrid %in% allsrrid) |>
  as.data.table() -> tbl_barcode_allsrrid


# src ---------------------------------------------------------------------

# function ----------------------------------------------------------------

fn_get_af_heatmap_celltype <- function(thevariant) {
  # return: a data.table contains af per celltype per srrid per gseid
  # thevariant <- "1004G>A"
  thepos <- gsub(
    pattern = "(\\d+)([A-Z])>([A-Z])",
    replacement = "\\1",
    x = thevariant
  ) |>
    as.character()
  theref <- gsub(
    pattern = "(\\d+)([A-Z])>([A-Z])",
    replacement = "\\2",
    x = thevariant
  )
  thealt <- gsub(
    pattern = "(\\d+)([A-Z])>([A-Z])",
    replacement = "\\3",
    x = thevariant
  )

  tbl_covall_allsrrid |>
    dplyr::select(gseid, srrid, barcode, base, all_of(thepos)) |>
    tidyr::pivot_wider(
      names_from = base,
      values_from = !!sym(thepos),
      values_fill = 0
    ) |>
    as.data.table() |>
    dplyr::left_join(
      tbl_barcode_allsrrid |>
        dplyr::select(gseid, srrid, barcode, celltype),
      by = c("gseid", "srrid", "barcode")
    ) -> covall_thevariant

  export(
    covall_thevariant,
    glue::glue(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/new-variant-cell/homo-hete/celllevel/{thevariant}.qs"
    )
  )

  covall_thevariant |>
    tidyr::nest(
      .by = c(gseid, srrid, celltype)
    ) |>
    dplyr::mutate(
      !!thevariant := purrr::map(
        .x = data,
        .f = \(.dt, thealt) {
          .dt[, .(A, C, G, T)] |>
            sum(na.rm = TRUE) -> total_depth
          if (total_depth == 0) {
            return(0)
          } else {
            af <- sum(.dt[, ..thealt], na.rm = TRUE) / total_depth
            return(af)
          }
        },
        thealt = thealt
      )
    ) |>
    dplyr::select(-data) -> af

  export(
    af,
    glue::glue(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/new-variant-cell/homo-hete/clusterlevel/{thevariant}.qs"
    )
  )
}

# body --------------------------------------------------------------------

# Create output directories if they don't exist
fs::dir_create(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/new-variant-cell/homo-hete/celllevel",
  recurse = TRUE
)
fs::dir_create(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/new-variant-cell/homo-hete/clusterlevel",
  recurse = TRUE
)

# Log processing start
logger::log_info("Processing {length(thevariants)} variants with 50 cores...")
logger::log_info("Variants: {paste(head(thevariants, 5), collapse = ', ')}...")

# Use parallel processing with 50 cores
results <- parallel::mclapply(
  X = thevariants,
  FUN = function(.variant) {
    tryCatch(
      {
        fn_get_af_heatmap_celltype(.variant)
        logger::log_success("Completed variant: {.variant}")
        return(TRUE)
      },
      error = function(e) {
        logger::log_error("Failed to process variant {.variant}: {e$message}")
        return(FALSE)
      }
    )
  },
  mc.cores = 50
)

# Check results
successful <- sum(unlist(results), na.rm = TRUE)
failed <- length(thevariants) - successful

logger::log_info(
  "Processing completed: {successful} successful, {failed} failed"
)

if (failed > 0) {
  failed_variants <- thevariants[!unlist(results)]
  logger::log_warn("Failed variants: {paste(failed_variants, collapse = ', ')}")
}

# Close database connections
DBI::dbDisconnect(conn, shutdown = TRUE)
DBI::dbDisconnect(conn_all_hetero_af, shutdown = TRUE)

logger::log_success("All processing completed and connections closed")

# footer ------------------------------------------------------------------

# save image --------------------------------------------------------------
