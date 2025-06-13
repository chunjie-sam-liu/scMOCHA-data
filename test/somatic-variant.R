thegseid <- "GSE171555"
thesrrid <- "GSM5227130"
thevariant <- "12501G>A"
thepos <- "12501"

all_variant_cell_table |>
  dplyr::filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  dplyr::collect() ->
seid_srrid_variant

filename <- file.path(
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/TABLES",
  "cell_cov.{thegseid}_{thesrrid}.csv" |> glue()
)


cell_cov <- data.table::fread(
  filename,
  header = TRUE,
)

cell_cov |>
  dplyr::select(
    base, barcode, dplyr::all_of(thepos)
  ) ->
cell_cov_sel

cell_cov_sel |>
  dplyr::group_by(
    barcode
  ) |>
  dplyr::summarise(
    depth = sum(`12501`)
  ) ->
cell_cov_depth


cell_cov_sel |>
  dplyr::filter(
    base %in% c("G", "A")
  ) |>
  dplyr::mutate(
    base = dplyr::recode(
      base,
      "G" = "ref",
      "A" = "alt"
    )
  ) |>
  tidyr::pivot_wider(
    names_from = base,
    values_from = `12501`
  ) |>
  dplyr::left_join(
    cell_cov_depth,
    by = "barcode"
  ) ->
cell_cov_variant

cell_cov_variant |>
  dplyr::mutate(
    variant_type_sanity = dplyr::case_when(
      depth == 0 ~ "white",
      depth < 10 ~ "grey",
      depth > 10 & alt < 4 ~ "black",
      depth > 10 & alt >= 4 ~ "colorful"
    )
  ) |>
  dplyr::left_join(
    seid_srrid_variant |> dplyr::select(
      barcode, variant_type
    ),
    by = c("barcode")
  ) |>
  dplyr::filter(
    variant_type_sanity != variant_type
  )
