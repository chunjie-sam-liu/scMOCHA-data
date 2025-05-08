gse_data <- readr::read_rds(
  "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_data.rds"
)


gse_data |>
  dplyr::select(gseid, srrid, depth_read) |>
  tidyr::unnest(cols = depth_read) |>
  dplyr::group_by(pos) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = T),
  ) ->
gse_data_coverage


gse_data_coverage |>
  data.table::fwrite(
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/", "gse_data_coverage.csv"
    ),
    row.names = F
  )
