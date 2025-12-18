pcc <- import(
  file = "https://raw.githubusercontent.com/chunjie-sam-liu/chunjie-sam-liu.life/master/public/data/pcc.tsv"
) |>
  dplyr::arrange(cancer_types)


METADATA <- import(
  path(
    "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data",
    "gse_dataset_metadata_full.qs"
  )
) |>
  dplyr::mutate(
    Haplogroup_s = purrr::map_chr(
      .x = Haplogroup,
      .f = \(.x) {
        # if (stringr::str_starts(.x, "L")) {
        #   gsub("L", "L0", .x)
        # }
        gsub("\\d+.*", "", .x)
      }
    )
  ) |>
  dplyr::mutate(Sex = SEXPRED)

METADATA |>
  dplyr::select(Haplogroup_s) |>
  dplyr::distinct() |>
  dplyr::filter(Haplogroup_s != "") |>
  dplyr::arrange(Haplogroup_s) |>
  dplyr::mutate(color = pcc$color[c(1:21)])
