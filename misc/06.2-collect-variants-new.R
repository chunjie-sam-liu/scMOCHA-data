d <- import(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_dataset_metadata_full.qs"
)

basedir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw"
d |>
  dplyr::select(gseid) |>
  dplyr::distinct() |>
  dplyr::mutate(
    anno = path(
      basedir,
      gseid,
      "out",
      "{gseid}.scmocha.out.qs" |> glue::glue()
    )
  ) |>
  dplyr::filter(!file.exists(anno))
