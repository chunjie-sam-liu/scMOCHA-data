conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)


mt_nd1_ <- c("3460G>A", "3697G>A", "3634A>G", "3380G>A")

dplyr::tbl(
  conn_all_hetero_af,
  "allvariants"
) |>
  dplyr::filter(variant %in% mt_nd1_)
