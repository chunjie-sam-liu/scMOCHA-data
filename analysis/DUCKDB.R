conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
packageVersion("duckdb")

# DBI::dbDisconnect(conn_all_hetero_af, shutdown = TRUE)
tables <- DBI::dbListTables(conn_all_hetero_af)

# tables[grepl("fisher", tables)] |> purrr::map(
#   .f =~{
#     DBI::dbRemoveTable(conn_all_hetero_af, .x)
#   }
# )
tables
dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_bulk"
)

dplyr::tbl(
  conn_all_hetero_af,
  "allvariants_cell"
)
dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_cell"
)

tables

# DBI::dbExecute(
#   conn_all_hetero_af,
#   "CREATE OR REPLACE TABLE allvariants_af_cell AS SELECT * FROM allvariants_cell"
# )

DBI::dbListTables(conn_all_hetero_af)


dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_bulk"
)

dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_af_cluster"
)


dplyr::tbl(
  conn_all_hetero_af,
  "allvariants_cell_covered"
) |>
  as.data.table() -> tbl_allvariants


tbl_allvariants |>
  dplyr::select(
    -c(
      af,
      depth,
      variant_type,
      variant_in_cell_cluster,
      pos,
      barcode
    )
  ) |>
  dplyr::rename(
    barcode = celltype
  ) |>
  tidyr::nest(
    .by = c(gseid, srrid, barcode, variant)
  ) -> tbl_allvariants_grouped

#' Calculate allele frequency for a given variant and data
#' @param .x : variant
#' @param .y : data
#' data : tibble with counts
#' data: columns are   AFO   ARE   CFO   CRE   GFO   GRE   TFO   TRE
#' @return allele frequency
#' @example fn_cal_af("10398A>G", data)
fn_cal_af <- function(.x, .y) {
  .refalt <- strsplit(
    x = gsub("\\d+", "", .x),
    split = ">"
  )[[1]]
  .altcount <- .y |>
    dplyr::select(
      dplyr::contains(.refalt[2])
    )
  .refcount <- .y |>
    dplyr::select(
      dplyr::contains(.refalt[1])
    )

  .af <- sum(.altcount, na.rm = TRUE) / sum(.y, na.rm = TRUE)
  return(.af)
}

tbl_allvariants_grouped |>
  # head(200) |>
  dplyr::mutate(
    af = parallel::mcmapply(
      FUN = fn_cal_af,
      .x = variant,
      .y = data,
      mc.cores = 20
    )
  ) |>
  dplyr::select(-data) -> tbl_allvariants_grouped_with_af

DBI::dbWriteTable(
  conn_all_hetero_af,
  "allvariants_af_cluster",
  tbl_allvariants_grouped_with_af,
  overwrite = TRUE,
  append = FALSE,
  temporary = FALSE
)
DBI::dbListTables(conn_all_hetero_af)
rm(tbl_allvariants_grouped)
gc()


tbl_allvariants |>
  dplyr::select(
    -c(
      af,
      depth,
      variant_type,
      variant_in_cell_cluster,
      pos,
      barcode,
      celltype
    )
  ) |>
  tidyr::nest(
    .by = c(gseid, srrid, variant)
  ) -> tbl_allvariants_bulk


tbl_allvariants_bulk |>
  dplyr::mutate(
    barcode = "bulk"
  ) |>
  dplyr::mutate(
    af = parallel::mcmapply(
      FUN = fn_cal_af,
      .x = variant,
      .y = data,
      mc.cores = 20
    )
  ) |>
  dplyr::select(-data) -> tbl_allvariants_bulk_with_af

DBI::dbWriteTable(
  conn_all_hetero_af,
  "allvariants_af_bulk",
  tbl_allvariants_bulk_with_af,
  overwrite = TRUE,
  append = FALSE,
  temporary = FALSE
)

DBI::dbListTables(conn_all_hetero_af)
dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_altdepth_cell"
)

dplyr::tbl(
  conn_all_hetero_af,
  "all_hetero_sumdepth_cell"
)


tbl_allvariants |>
  dplyr::select(
    -c(
      af,
      depth,
      variant_type,
      variant_in_cell_cluster,
      pos
    )
  )
