gse_data <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_data.qs"
)


gse_data |>
  dplyr::select(gseid, srrid, chemistry, depth_read) |>
  tidyr::unnest(cols = depth_read) ->
gse_data_depth

gse_data_depth |>
  dplyr::group_by(pos) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = T),
  ) ->
gse_data_coverage


gse_data_coverage |>
  export(
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "gse_data_coverage.csv"
    ),
    format = "both"
  )


gse_data_depth |>
  dplyr::group_by(chemistry, pos) |>
  dplyr::summarise(
    depth = mean(depth, na.rm = T),
  ) |>
  dplyr::ungroup() ->
gse_data_coverage_chemistry

gse_data_coverage_chemistry |>
  export(
    file = file.path(
      "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/", "gse_data_coverage_chemistry.csv"
    ),
    format = "both"
  )
