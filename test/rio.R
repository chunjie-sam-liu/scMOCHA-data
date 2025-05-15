mtcars |> head()

mtcars |>
  head() |>
  dplyr::select(cyl, gear) |>
  dplyr::arrange(cyl) |>
  tidyr::expand(cyl, gear)

mtcars |>
  dplyr::select(cyl) |>
  dplyr::lag(cyl)
