allvariants <- import(
  "/home/cliu68/github/scMOCHA-data/high-res-MANUSCRIPTFIGURES/SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF-CLEAN.xlsx"
)


conflicted::conflict_prefer("filter", "dplyr")


allvariants |>
  filter(
    `Variant type` == "homoplasmic"
  ) |>
  filter(
    `Sample ID` == "GSM4509019"
  ) |>
  # filter(
  #   B < 0.95
  # ) |>
  filter(
    Variant == "9540T>C"
  ) |>
  View()


allvariants |>
  filter(
    `Variant type` == "homoplasmic"
  ) |>
  filter(
    `Sample ID` == "GSM4509019"
  ) |>
  filter(
    `CD4_T` < 0.95
  ) |>
  filter(
    Variant == "8414C>T"
  ) |>
  View()
