thevariant <- "3727T>C"
thegseid <- "GSE235050"
thesrrid <- "GSM7493833"

tbl_all_hetero_af_cell <- dplyr::tbl(
  conn,
  "all_hetero_af_cell"
)

tbl_allvariants_cell <- dplyr::tbl(
  conn,
  "allvariants_cell"
)


tbl_all_hetero_af_cell |>
  filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  select(
    -variant
  ) |>
  rename(
    af1 = af
  )

tbl_allvariants_cell |>
  filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  select(-variant_in_cell_cluster) |>
  select(-variant) |>
  inner_join(
    tbl_all_hetero_af_cell |>
      filter(
        gseid == thegseid,
        srrid == thesrrid,
        variant == thevariant
      ) |>
      select(
        -variant
      ) |>
      rename(
        af_old = af
      )
  ) |>
  as.data.table() |>
  filter(af != af_old)

m <- import(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/final/GSM7493833/variant_info_from_heatmap.qs"
)
m

DBI::dbListTables(conn)
tbl(conn, "all_hetero_af_cell")

m |> filter(variant == thevariant) |> as.data.table()

m |>
  filter(variant == thevariant) |>
  # filter(celltype == "B") |>
  ggplot(aes(x = af, fill = variant_type)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 30
  ) +
  facet_wrap(~celltype)


tbl_all_hetero_af_cell |>
  filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  select(
    -variant
  ) |>
  rename(
    af1 = af
  ) |>
  ggplot(aes(x = af1)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 30
  ) +
  facet_wrap(~celltype)


tbl_all_hetero_af_cell |>
  filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  select(
    -variant
  ) |>
  rename(
    af1 = af
  ) |>
  filter(af1 == 1) |>
  as.data.table() -> allcells1


tbl_allvariants_cell |>
  filter(
    gseid == thegseid,
    srrid == thesrrid,
    variant == thevariant
  ) |>
  select(-variant_in_cell_cluster) |>
  select(-variant) |>
  as.data.table() |>
  slice(
    match(allcells1$barcode, barcode)
  )


m |>
  filter(
    variant == thevariant
  ) |>
  slice(
    match(allcells1$barcode, barcode)
  ) |>
  count(variant_type)
allcells1
