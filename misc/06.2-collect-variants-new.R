d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)


d |>
  dplyr::select(gseid, srrdir) |>
  dplyr::mutate(gsedir = dirname(dirname(dirname(srrdir)))) |>
  dplyr::select(gseid, gsedir) |>
  dplyr::distinct() |>
  dplyr::mutate(
    newrun = path(
      gsedir,
      gseid,
      "final",
      "out",
      "{gseid}.scmocha.out.qs" |> glue::glue()
    )
  )
